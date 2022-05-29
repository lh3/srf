#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <assert.h>
#include "ketopt.h"
#include "rld0.h"

#define Calloc(type, cnt) ((type*)calloc((cnt), sizeof(type)))

typedef struct {
	int64_t k, l;
	int d, c;
} pair64_t;

typedef struct {
	rld_t *e;
	uint64_t s_top;
	pair64_t *stack;
	uint64_t ok[6], ol[6];
	pair64_t top;
} aux_t;

int main(int argc, char *argv[])
{
	int c, min_occ = 100, depth = 51, i, n = 0, use_mmap = 0;
	pair64_t *p;
	aux_t *aux;
	char *str;
	ketopt_t o = KETOPT_INIT;

	while ((c = ketopt(&o, argc, argv, 1, "k:m:M", 0)) >= 0) {
		if (c == 'k') depth = atol(o.arg);
		else if (c == 'm') min_occ = atol(o.arg);
		else if (c == 'M') use_mmap = 1;
	}
	if (o.ind == argc) {
		fprintf(stderr, "Usage: fmd-occ [options] <in1.fmd> [in2.fmd [...]]\n");
		fprintf(stderr, "Options:\n");
		fprintf(stderr, "  -k INT       k-mer length [%d]\n", depth);
		fprintf(stderr, "  -m INT       min k-mer occurrence [%d]\n", min_occ);
		return 1;
	}

	n = argc - o.ind;
	aux = Calloc(aux_t, n);
	for (i = 0; i < n; ++i) {
		aux_t *q = &aux[i];
		q->e = use_mmap? rld_restore_mmap(argv[o.ind + i]) : rld_restore(argv[o.ind + i]);
		assert(q->e);
		q->stack = Calloc(pair64_t, depth + 1);
		p = &q->stack[q->s_top++];
		p->k = 0, p->l = q->e->mcnt[0], p->d = p->c = 0;
	}
	str = Calloc(char, depth + 1);
	str[depth] = 0;
	while (1) {
		int a;
		for (i = 0; i < n; ++i) {
			if (aux[i].s_top == 0) break;
			aux[i].top = aux[i].stack[--aux[i].s_top];
		}
		if (i < n) break;
		if (aux->top.d > 0) str[depth - aux->top.d] = "$ACGTN"[aux->top.c];
		for (i = 0; i < n; ++i)
			rld_rank2a(aux[i].e, aux[i].top.k, aux[i].top.l, aux[i].ok, aux[i].ol);
		for (a = 1; a <= 4; ++a) {
			for (i = 0; i < n; ++i)
				if (aux[i].ol[a] - aux[i].ok[a] >= min_occ) break;
			if (i == n) continue;
			str[depth - aux->top.d - 1] = "$ACGTN"[a];
			if (aux->top.d != depth - 1) {
				for (i = 0; i < n; ++i) {
					aux_t *q = &aux[i];
					p = &q->stack[q->s_top++];
					p->k = q->e->cnt[a] + q->ok[a];
					p->l = q->e->cnt[a] + q->ol[a];
					p->d = aux->top.d + 1;
					p->c = a;
				}
			} else {
				fputs(str, stdout);
				for (i = 0; i < n; ++i)
					printf("\t%ld", (long)(aux[i].ol[a] - aux[i].ok[a]));
				putchar('\n');
			}
		}
	}
	free(str);
	for (i = 0; i < n; ++i) {
		free(aux[i].stack);
		rld_destroy(aux[i].e);
	}
	return 0;
}
