#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <assert.h>
#include <stdio.h>
#include <zlib.h>
#include "khashl.h"
#include "ketopt.h"
#include "kseq.h"
#include "ksort.h"

#define Malloc(type, cnt)       ((type*)malloc((cnt) * sizeof(type)))
#define Calloc(type, cnt)       ((type*)calloc((cnt), sizeof(type)))
#define Realloc(type, ptr, cnt) ((type*)realloc((ptr), (cnt) * sizeof(type)))

KSTREAM_INIT(gzFile, gzread, 0x10000)

typedef struct {
	uint32_t cnt;
	int32_t len;
	uint32_t k, used;
	uint8_t *seq[2];
} ca_kmer_t;

static inline int32_t ca_kmer_eq(const ca_kmer_t x, const ca_kmer_t y)
{
	return (x.len == y.len && memcmp(x.seq[0], y.seq[0], x.len) == 0);
}

static inline uint32_t ca_kmer_hash(const ca_kmer_t x)
{
	return kh_hash_bytes(x.len, x.seq[0]);
}

KHASHL_CSET_INIT(, ca_kh_t, ca_kh, ca_kmer_t, ca_kmer_hash, ca_kmer_eq)

#define kmer_key(x) ((x).cnt)
KRADIX_SORT_INIT(ca_kmer, ca_kmer_t, kmer_key, 4)

unsigned char seq_nt4_table[256] = {
	0, 1, 2, 3,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  3, 3, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  3, 3, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};

ca_kh_t *ca_kmer_read(const char *fn)
{
	ca_kh_t *h = 0;
	gzFile fp;
	kstream_t *ks;
	kstring_t str = {0,0,0};
	int32_t dret;
	ca_kmer_t kmer;
	uint8_t *swap;

	fp = gzopen(fn, "rb");
	if (fp == 0) return 0;
	ks = ks_init(fp);
	h = ca_kh_init();

	kmer.len = -1, kmer.seq[0] = kmer.seq[1] = 0;
	while (ks_getuntil(ks, KS_SEP_LINE, &str, &dret) >= 0) {
		int32_t i, absent;
		long cnt;
		char *p;
		khint_t k;
		for (i = 0; i < str.l; ++i)
			if (str.s[i] == '\t' || str.s[i] == ' ')
				break;
		assert(i + 1 < str.l);
		cnt = strtol(&str.s[i+1], &p, 10);
		if (cnt <= 0) continue;
		kmer.cnt = cnt;
		if (kmer.len < 0) {
			assert((i&1) == 1);
			kmer.len = i;
			kmer.seq[0] = Malloc(uint8_t, kmer.len * 2);
			kmer.seq[1] = kmer.seq[0] + kmer.len;
			swap = Malloc(uint8_t, kmer.len);
		}
		assert(kmer.len == i);
		for (i = 0; i < kmer.len; ++i) {
			int32_t c = seq_nt4_table[(uint8_t)str.s[i]];
			if (c >= 4) break;
			kmer.seq[0][i] = c;
			kmer.seq[1][kmer.len - 1 - i] = 3 - c;
		}
		if (i < kmer.len) continue;
		if (memcmp(kmer.seq[0], kmer.seq[1], kmer.len) > 0) {
			memcpy(swap, kmer.seq[0], kmer.len);
			memcpy(kmer.seq[0], kmer.seq[1], kmer.len);
			memcpy(kmer.seq[1], swap, kmer.len);
		}
		k = ca_kh_put(h, kmer, &absent);
		if (absent) {
			ca_kmer_t *q = &kh_key(h, k);
			q->seq[0] = Malloc(uint8_t, kmer.len * 2);
			q->seq[1] = q->seq[0] + kmer.len;
			memcpy(q->seq[0], kmer.seq[0], kmer.len * 2);
		}
	}
	free(str.s);
	free(swap);
	free(kmer.seq[0]);
	ks_destroy(ks);
	gzclose(fp);
	return h;
}

void ca_gen(ca_kh_t *h)
{
	int32_t i, n;
	khint_t k;
	ca_kmer_t *a;
	n = kh_size(h);
	a = Malloc(ca_kmer_t, n);
	for (k = 0, i = 0; k != kh_end(h); ++k)
		if (kh_exist(h, k))
			kh_key(h, k).k = k, kh_key(h, k).used = 0, a[i++] = kh_key(h, k);
	radix_sort_ca_kmer(a, a + n);
	for (i = 0; i < n>>1; ++i) {
		ca_kmer_t t = a[i];
		a[i] = a[n - 1 - i], a[n - 1 - i] = t;
	}
	for (i = 0; i < n; ++i) {
		printf("%d\n", a[i].cnt);
	}
	free(a);
}

int main(int argc, char *argv[])
{
	int32_t c;
	ketopt_t o = KETOPT_INIT;
	ca_kh_t *h;
	while ((c = ketopt(&o, argc, argv, 1, "", 0)) >= 0) {
	}
	if (o.ind == argc) {
		fprintf(stderr, "Usage: circasm <in.txt>\n");
		return 1;
	}
	h = ca_kmer_read(argv[o.ind]);
	fprintf(stderr, "[M::%s] read %d distinct k-mers\n", __func__, kh_size(h));
	ca_gen(h);
	return 0;
}
