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

#define generic_key(x) (x)
KRADIX_SORT_INIT(ca64, uint64_t, generic_key, 8)

typedef struct {
	uint32_t cnt;
	int32_t len;
	uint32_t k, flag;
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

void ca_kmer_canonical(ca_kmer_t *t, uint8_t *swap)
{
	if (memcmp(t->seq[0], t->seq[1], t->len) > 0) {
		memcpy(swap, t->seq[0], t->len);
		memcpy(t->seq[0], t->seq[1], t->len);
		memcpy(t->seq[1], swap, t->len);
	}
}

void ca_kmer_append(ca_kmer_t *t, int32_t len, const uint8_t *seq, int32_t c, uint8_t *swap)
{
	int32_t i;
	t->len = len;
	memcpy(t->seq[0], seq, len - 1);
	t->seq[0][len - 1] = c;
	for (i = 0; i < len; ++i)
		t->seq[1][len - i - 1] = 3 - t->seq[0][i];
	ca_kmer_canonical(t, swap);
}

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
	uint8_t *swap = 0;

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
		ca_kmer_canonical(&kmer, swap);
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
	int32_t i, n, l_seq, m_seq;
	uint64_t *a;
	khint_t k;
	ca_kmer_t t;
	uint8_t *seq, *swap;

	n = kh_size(h);
	if (n == 0) return;
	a = Malloc(uint64_t, n);
	for (k = 0, i = 0; k != kh_end(h); ++k)
		if (kh_exist(h, k))
			kh_key(h, k).k = k, kh_key(h, k).flag = (uint32_t)-1, a[i++] = (uint64_t)kh_key(h, k).cnt<<32 | k;
	radix_sort_ca64(a, a + n);
	for (i = 0; i < n>>1; ++i) { // change to the descending order
		uint64_t t = a[i];
		a[i] = a[n - 1 - i], a[n - 1 - i] = t;
	}

	t = kh_key(h, (uint32_t)a[0]);
	t.seq[0] = Calloc(uint8_t, t.len * 2);
	t.seq[1] = t.seq[0] + t.len;
	swap = Calloc(uint8_t, t.len);

	l_seq = m_seq = 0;
	seq = 0;
	for (i = 0; i < n; ++i) {
		khint_t k0 = (uint32_t)a[i];
		ca_kmer_t *q = &kh_key(h, k0);
		int32_t done;

		if (q->flag != (uint32_t)-1) continue;
		q->flag = k0;

		if (q->len >= m_seq) {
			m_seq = q->len * 2;
			seq = Calloc(uint8_t, m_seq);
		}
		memcpy(seq, q->seq[0], q->len);
		l_seq = q->len;
		done = 0;
		while (1) {
			int32_t cnt[4], max, c, max_c;
			khint_t kk[4];
			for (c = 0; c < 4; ++c) {
				khint_t k;
				ca_kmer_append(&t, q->len, &seq[l_seq - q->len + 1], c, swap);
				kk[c] = k = ca_kh_get(h, t);
				cnt[c] = k != kh_end(h)? kh_key(h, k).cnt : 0;
			}
			for (c = 0, max = 0, max_c = -1; c < 4; ++c)
				if (cnt[c] > max)
					max = cnt[c], max_c = c;
			if (max == 0) break;
			ca_kmer_append(&t, q->len, &seq[l_seq - q->len + 1], max_c, swap);
			if (ca_kmer_eq(kh_key(h, k0), t)) {
				done = 1;
				break;
			}
			for (c = 0; c < 4; ++c)
				if (kk[c] != kh_end(h))
					kh_key(h, kk[c]).flag = k0;
			if (l_seq == m_seq) {
				m_seq += (m_seq>>1) + 16;
				seq = Realloc(uint8_t, seq, m_seq);
			}
			seq[l_seq++] = max_c;
			q = &kh_key(h, kk[max_c]);
		}
		if (done) {
			int32_t j;
			for (j = 0; j < l_seq; ++j)
				seq[j] = "ACGT"[seq[j]];
			fwrite(seq, 1, l_seq, stdout);
			putchar('\n');
		}
	}
	free(seq);
	free(a);
	free(swap);
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
