/* substitute_filter.c  -- fast inner engine for substitute_fragments.py
 *
 * stdin  : decompressed fragment TSV (chrom\tstart\tstop\t<16mer>-<suf>\tcount\n)
 * stdout : filtered stream (affected barcodes dropped; substitute lines get relabeled
 *          copies appended immediately after, preserving sort order)
 * stderr : periodic progress + final read/written counts
 * argv[1]: map file with lines:
 *            A\t<16mer>                  affected barcode (drop its fragments)
 *            S\t<sub16mer>\t<aff16mer>   substitute -> affected mapping
 *
 * 16-mers are pure ACGT and exactly 16 long, so each encodes to a uint32 (2 bits/base).
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <stdint.h>

static int8_t B2[256];

static inline int enc16(const char *p, uint32_t *out) {
    uint32_t k = 0;
    for (int i = 0; i < 16; i++) {
        int8_t v = B2[(unsigned char)p[i]];
        if (v < 0) return 0;
        k = (k << 2) | (uint32_t)v;
    }
    *out = k;
    return 1;
}
static inline void dec16(uint32_t k, char *p) {
    static const char L[4] = {'A', 'C', 'G', 'T'};
    for (int i = 15; i >= 0; i--) { p[i] = L[k & 3]; k >>= 2; }
}

#define HSZ  4096
#define HMASK (HSZ - 1)
static uint32_t aff_key[HSZ];  static uint8_t aff_used[HSZ];
static uint32_t sub_key[HSZ];  static uint8_t sub_used[HSZ];
static uint32_t *sub_list[HSZ]; static int sub_n[HSZ], sub_cap[HSZ];

static void aff_add(uint32_t k) {
    uint32_t h = (k * 2654435761u) & HMASK;
    while (aff_used[h]) { if (aff_key[h] == k) return; h = (h + 1) & HMASK; }
    aff_used[h] = 1; aff_key[h] = k;
}
static inline int aff_has(uint32_t k) {
    uint32_t h = (k * 2654435761u) & HMASK;
    while (aff_used[h]) { if (aff_key[h] == k) return 1; h = (h + 1) & HMASK; }
    return 0;
}
static void sub_add(uint32_t s, uint32_t a) {
    uint32_t h = (s * 2654435761u) & HMASK;
    while (sub_used[h]) { if (sub_key[h] == s) break; h = (h + 1) & HMASK; }
    if (!sub_used[h]) { sub_used[h] = 1; sub_key[h] = s; sub_cap[h] = 4; sub_n[h] = 0;
                        sub_list[h] = malloc(4 * sizeof(uint32_t)); }
    if (sub_n[h] == sub_cap[h]) { sub_cap[h] *= 2;
                        sub_list[h] = realloc(sub_list[h], sub_cap[h] * sizeof(uint32_t)); }
    sub_list[h][sub_n[h]++] = a;
}
static inline int sub_find(uint32_t s) {
    uint32_t h = (s * 2654435761u) & HMASK;
    while (sub_used[h]) { if (sub_key[h] == s) return (int)h; h = (h + 1) & HMASK; }
    return -1;
}

static void write_all(const char *p, size_t n) {
    while (n) { ssize_t w = write(1, p, n); if (w <= 0) { perror("write"); exit(1); } p += w; n -= w; }
}

static long long g_nr = 0, g_nw = 0;   /* lines read / written */

/* Process one line [p, nl). Appends to out (*olen). Updates counters. */
static inline void proc(const char *p, const char *nl, char *out, size_t *olen) {
    size_t len = (size_t)(nl - p);
    const char *bc = NULL; int tabs = 0;
    for (const char *c = p; c < nl; c++)
        if (*c == '\t') { if (++tabs == 3) { bc = c + 1; break; } }

    g_nr++;
    uint32_t k;
    if (bc && bc + 16 <= nl && enc16(bc, &k)) {
        if (aff_has(k)) return;                       /* affected: drop */
        memcpy(out + *olen, p, len); *olen += len; out[(*olen)++] = '\n'; g_nw++;
        int hi = sub_find(k);
        if (hi >= 0) {                                /* substitute: append copies */
            size_t prelen = (size_t)(bc - p);
            const char *rest = bc + 16;
            size_t restlen = (size_t)(nl - rest);
            for (int i = 0; i < sub_n[hi]; i++) {
                memcpy(out + *olen, p, prelen); *olen += prelen;
                dec16(sub_list[hi][i], out + *olen); *olen += 16;
                memcpy(out + *olen, rest, restlen); *olen += restlen;
                out[(*olen)++] = '\n'; g_nw++;
            }
        }
    } else {                                          /* unparseable: passthrough */
        memcpy(out + *olen, p, len); *olen += len; out[(*olen)++] = '\n'; g_nw++;
    }
}

int main(int argc, char **argv) {
    if (argc < 2) { fprintf(stderr, "usage: %s mapfile < in > out\n", argv[0]); return 2; }
    memset(B2, -1, sizeof(B2));
    B2['A'] = 0; B2['C'] = 1; B2['G'] = 2; B2['T'] = 3;

    FILE *mf = fopen(argv[1], "r");
    if (!mf) { perror("mapfile"); return 1; }
    char *ln = NULL; size_t lc = 0; ssize_t lr;
    while ((lr = getline(&ln, &lc, mf)) > 0) {
        if (ln[0] == 'A' && ln[1] == '\t') {
            uint32_t k; if (enc16(ln + 2, &k)) aff_add(k);
        } else if (ln[0] == 'S' && ln[1] == '\t') {
            uint32_t s, a;
            if (enc16(ln + 2, &s) && enc16(ln + 2 + 16 + 1, &a)) sub_add(s, a);
        }
    }
    free(ln); fclose(mf);

    const size_t CAP  = (size_t)1 << 26;    /* 64 MB input buffer  */
    const size_t OCAP = (size_t)1 << 26;    /* 64 MB output buffer */
    char *buf = malloc(CAP);
    char *out = malloc(OCAP);
    size_t have = 0, olen = 0;
    long long next_report = 50000000LL;

    for (;;) {
        ssize_t r = read(0, buf + have, CAP - have);
        if (r < 0) { perror("read"); return 1; }
        have += (size_t)r;

        char *p = buf, *end = buf + have, *nl;
        while ((nl = memchr(p, '\n', (size_t)(end - p)))) {
            proc(p, nl, out, &olen);
            if (olen > OCAP - (1 << 16)) { write_all(out, olen); olen = 0; }
            if (g_nr >= next_report) {
                fprintf(stderr, "\r  %lldM read  %lldM written", g_nr / 1000000, g_nw / 1000000);
                next_report += 50000000LL;
            }
            p = nl + 1;
        }
        have = (size_t)(end - p);
        memmove(buf, p, have);
        if (r == 0) break;
        if (have == CAP) { fprintf(stderr, "\nline exceeds buffer\n"); return 1; }
    }
    if (have > 0) {                                   /* trailing line, no newline */
        proc(buf, buf + have, out, &olen);
    }
    if (olen) write_all(out, olen);
    fprintf(stderr, "\r  %lld read, %lld written\n", g_nr, g_nw);
    return 0;
}
