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

static inline uint32_t murmur_32_scramble(uint32_t k) {
    k *= 0xcc9e2d51;
    k = (k << 15) | (k >> 17);
    k *= 0x1b873593;
    return k;
}

uint32_t murmur3_32(const uint8_t* key, size_t len, uint32_t seed)
{
	uint32_t h = seed;
    uint32_t k;
	size_t i;
    for (i = len >> 2; i; i--) {
        memcpy(&k, key, sizeof(uint32_t));
        key += sizeof(uint32_t);
        h ^= murmur_32_scramble(k);
        h = (h << 13) | (h >> 19);
        h = h * 5 + 0xe6546b64;
    }
    k = 0;
    for (i = len & 3; i; i--) {
        k <<= 8;
        k |= key[i - 1];
    }
    h ^= murmur_32_scramble(k);
	h ^= len;
	h ^= h >> 16;
	h *= 0x85ebca6b;
	h ^= h >> 13;
	h *= 0xc2b2ae35;
	h ^= h >> 16;
	return h;
}

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
	return murmur3_32(x.seq[0], x.len, 11);
}

KHASHL_CSET_INIT(, ca_kh_t, ca_kh, ca_kmer_t, ca_kmer_hash, ca_kmer_eq)

int32_t ca_kmer_canonical(ca_kmer_t *t, uint8_t *swap)
{
	if (memcmp(t->seq[0], t->seq[1], t->len) > 0) {
		memcpy(swap, t->seq[0], t->len);
		memcpy(t->seq[0], t->seq[1], t->len);
		memcpy(t->seq[1], swap, t->len);
		return 1;
	} else return 0;
}

int32_t ca_kmer_append(ca_kmer_t *t, int32_t len, const uint8_t *seq, int32_t c, uint8_t *swap)
{
	int32_t i;
	t->len = len;
	memcpy(t->seq[0], seq, len - 1);
	t->seq[0][len - 1] = c;
	for (i = 0; i < len; ++i)
		t->seq[1][len - i - 1] = 3 - t->seq[0][i];
	return ca_kmer_canonical(t, swap);
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
		} else kh_key(h, k).cnt += cnt;
	}
	free(str.s);
	free(swap);
	free(kmer.seq[0]);
	ks_destroy(ks);
	gzclose(fp);
	return h;
}

void ca_gen(ca_kh_t *h, const char *prefix)
{
	int32_t i, n, l_seq, m_seq, n_circ = 0;
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
		int32_t done, min_cnt, max_cnt;

		if (q->flag != (uint32_t)-1) continue;

		if (q->len >= m_seq) {
			m_seq = q->len * 2;
			seq = Calloc(uint8_t, m_seq);
		}
		memcpy(seq, q->seq[0], q->len);
		l_seq = q->len;
		min_cnt = max_cnt = q->cnt;
		done = 0;
		while (1) {
			int32_t cnt[4], max, c, max_c;
			khint_t kk[4];
			q->flag = k0;
			min_cnt = min_cnt < q->cnt? min_cnt : q->cnt;
			max_cnt = max_cnt > q->cnt? max_cnt : q->cnt;
			//fprintf(stderr, "X\t%d\t%d\t%d\n", i, q->cnt, l_seq);
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
			if (l_seq == m_seq) {
				m_seq += (m_seq>>1) + 16;
				seq = Realloc(uint8_t, seq, m_seq);
			}
			seq[l_seq++] = max_c;
			q = &kh_key(h, kk[max_c]);
			if (q->flag != (uint32_t)-1) break;
		}
		if (done) {
			int32_t j;
			++n_circ;
			putchar('>');
			if (prefix) printf("%s#", prefix);
			printf("circ%d-%d min=%d,max=%d\n", n_circ, l_seq - q->len + 1, min_cnt, max_cnt);
			for (j = 0; j < l_seq - q->len + 1; ++j)
				seq[j] = "ACGT"[seq[j]];
			fwrite(seq, 1, l_seq - q->len + 1, stdout);
			putchar('\n');
		}
	}
	free(seq);
	free(a);
	free(swap);
}

typedef struct {
	uint32_t par, k0;
	int32_t w;
} ninfo_t;

typedef struct {
	uint32_t k;
	int32_t cnt;
	int32_t w;
} elem_t;

typedef struct {
	int32_t n, m;
	elem_t *a;
} heap_t;

#define CA_PAR_UNSET ((uint32_t)-1)
#define CA_PAR_START ((uint32_t)-2)

#define elem_lt(x, y) ((x).cnt < (y).cnt)
KSORT_INIT(ca_elem, elem_t, elem_lt)

static inline void heap_insert(heap_t *hp, const ca_kh_t *h, uint32_t k, int32_t w)
{
	elem_t *p;
	if (hp->n == hp->m) {
		hp->m += (hp->m>>1) + 16;
		hp->a = Realloc(elem_t, hp->a, hp->m);
	}
	p = &hp->a[hp->n++];
	p->k = k, p->w = w, p->cnt = kh_key(h, k).cnt;
	ks_heapup_ca_elem(hp->n, hp->a);
}

static inline elem_t heap_extract_max(heap_t *hp)
{
	elem_t t;
	t = hp->a[0]; // FIXME: not working if p->n == 0
	hp->a[0] = hp->a[--(hp->n)];
	ks_heapdown_ca_elem(0, hp->n, hp->a);
	return t;
}

void ca_gen_heap(const ca_kh_t *h, const char *prefix)
{
	int32_t i, n, len, m_seq = 0;
	uint64_t *a;
	ninfo_t *f;
	khint_t k;
	uint8_t *swap;
	char *seq = 0;
	ca_kmer_t tmp;
	heap_t hp = {0,0,0};

	n = kh_size(h);
	if (n == 0) return;

	// collect count array a[]
	a = Malloc(uint64_t, n);
	for (k = 0, i = 0; k != kh_end(h); ++k)
		if (kh_exist(h, k))
			kh_key(h, k).k = k, kh_key(h, k).flag = (uint32_t)-1, a[i++] = (uint64_t)kh_key(h, k).cnt<<32 | k;

	// sort by the descending order
	radix_sort_ca64(a, a + n);
	for (i = 0; i < n>>1; ++i) { // change to the descending order
		uint64_t t = a[i];
		a[i] = a[n - 1 - i], a[n - 1 - i] = t;
	}

	tmp = kh_key(h, (uint32_t)a[0]);
	len = tmp.len;
	tmp.seq[0] = Calloc(uint8_t, len * 2);
	tmp.seq[1] = tmp.seq[0] + len;
	swap = Calloc(uint8_t, len);
	f = Calloc(ninfo_t, kh_end(h));
	for (i = 0; i < kh_end(h); ++i)
		f[i].par = CA_PAR_UNSET, f[i].w = -1, f[i].k0 = CA_PAR_UNSET;

	for (i = 0; i < n; ++i) {
		khint_t k0 = (uint32_t)a[i];
		int32_t succ = 0;
		if (f[k0].par != CA_PAR_UNSET) continue;
		hp.n = 0;
		heap_insert(&hp, h, k0, 0);
		f[k0].par = CA_PAR_START, f[k0].w = 0;
		while (hp.n > 0) {
			elem_t e;
			ca_kmer_t *q;
			int32_t c;
			e = heap_extract_max(&hp);
			if (e.k == k0 && f[k0].par != CA_PAR_START) {
				succ = 1;
				break;
			}
			q = &kh_key(h, e.k);
			for (c = 0; c < 4; ++c) {
				int32_t w = ca_kmer_append(&tmp, q->len, &q->seq[e.w][1], c, swap);
				khint_t k = ca_kh_get(h, tmp);
				if (k == kh_end(h)) continue;
				if (f[k].par == CA_PAR_UNSET || (f[k].par == CA_PAR_START && f[k].w == 0)) {
					//if (k0 == 1137185) printf("X\t%d\t%d:%d -> %d:%d\n", k0, e.k, e.w, k, w);
					heap_insert(&hp, h, k, w);
					f[k].par = e.k, f[k].w = e.w, f[k].k0 = k0;
				}
			}
		}
		if (succ) {
			khint_t k = f[k0].par;
			int32_t l = 0, w = f[k0].w, t;
			//printf("==> start: %d\n", k0);
			while (1) {
				//printf("%d\t%d\n", k, f[k].k0);
				if (l == m_seq) {
					m_seq += (m_seq>>1) + 16;
					seq = Realloc(char, seq, m_seq);
				}
				seq[l++] = "ACGT"[kh_key(h, k).seq[w][0]];
				if (k == k0) break;
				w = f[k].w, k = f[k].par;
			}
			seq[l] = 0;
			for (i = 0; i < l>>1; ++i)
				t = seq[i], seq[i] = seq[l-i-1], seq[l-i-1] = t;
			puts(seq);
		}
	}
	free(f);
	free(swap);
	free(tmp.seq[0]);
	free(a);
}

int main(int argc, char *argv[])
{
	int32_t c;
	ketopt_t o = KETOPT_INIT;
	ca_kh_t *h;
	char *prefix = 0;
	while ((c = ketopt(&o, argc, argv, 1, "p:", 0)) >= 0) {
		if (c == 'p') prefix = o.arg;
	}
	if (o.ind == argc) {
		fprintf(stderr, "Usage: srf [options] <in.txt>\n");
		fprintf(stderr, "Options:\n");
		fprintf(stderr, "  -p STR     output prefix []\n");
		return 1;
	}
	h = ca_kmer_read(argv[o.ind]);
	fprintf(stderr, "[M::%s] read %d distinct k-mers\n", __func__, kh_size(h));
#if 1
	ca_gen_heap(h, prefix);
#else
	ca_gen(h, prefix);
#endif
	return 0;
}
