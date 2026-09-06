/* Optional experimental backend, built against OpenZL v0.2.0. Not linked into aie.
 * Usage: pathbench-openzl MODE LEVEL INPUT OUTPUT
 * MODE 0 numeric generic, 1 delta, 2 field-LZ, 3 byte zstd, 4 range/generic,
 * 5 range/FSE, 6 range/zstd, 7 flatpack, 8 tokenize/FSE, 9 bitpack.
 * Every persisted frame is reopened and verified. Times exclude file I/O.
 */
#define _POSIX_C_SOURCE 200809L
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "openzl/zl_compressor.h"
#include "openzl/zl_decompress.h"
#include "openzl/zl_version.h"
#include "zstd.h"
static int mode, level, empty;
static double now(void) { struct timespec t; clock_gettime(CLOCK_MONOTONIC,&t); return t.tv_sec+t.tv_nsec*1e-9; }
static void fail(const char *s) {fprintf(stderr,"%s\n",s);exit(1);}
static size_t checked(ZL_Report r) {if(ZL_isError(r)) fail("OpenZL codec error");return ZL_validResult(r);}
static ZL_GraphID graph(ZL_Compressor *c) {
    checked(ZL_Compressor_setParameter(c,ZL_CParam_formatVersion,ZL_MAX_FORMAT_VERSION));
    checked(ZL_Compressor_setParameter(c,ZL_CParam_compressionLevel,level));
    ZL_GraphID g=mode==2?ZL_GRAPH_FIELD_LZ:ZL_GRAPH_COMPRESS_GENERIC;
    if(mode==1) g=ZL_Compressor_registerStaticGraph_fromNode1o(c,ZL_NODE_DELTA_INT,g);
    if(mode==4) g=ZL_Compressor_registerStaticGraph_fromNode1o(c,ZL_NODE_RANGE_PACK,g);
    if(mode==5 || mode==6) {
        g=ZL_Compressor_registerStaticGraph_fromNode1o(c,ZL_NODE_CONVERT_NUM_TO_SERIAL_LE,mode==5?ZL_GRAPH_FSE:ZL_Compressor_registerZstdGraph_withLevel(c,level));
        g=ZL_Compressor_registerStaticGraph_fromNode1o(c,ZL_NODE_RANGE_PACK,g);
    }
    if(mode==7) g=ZL_GRAPH_FLATPACK;
    if(mode==8) {
        ZL_GraphID indices=ZL_Compressor_registerStaticGraph_fromNode1o(c,ZL_NODE_CONVERT_NUM_TO_SERIAL_LE,ZL_GRAPH_FSE);
        indices=ZL_Compressor_registerStaticGraph_fromNode1o(c,ZL_NODE_RANGE_PACK,indices);
        g=ZL_Compressor_registerTokenizeGraph(c,ZL_Type_numeric,0,ZL_GRAPH_BITPACK,indices);
    }
    if(mode==9) g=ZL_GRAPH_BITPACK;
    if(empty) g=ZL_GRAPH_STORE;
    return ZL_Compressor_registerStaticGraph_fromNode1o(c,ZL_NODE_INTERPRET_AS_LE64,g);
}
static unsigned char *read_file(const char *name,size_t *n) {
    FILE *f=fopen(name,"rb");if(!f)fail("cannot open input");
    if(fseek(f,0,SEEK_END))fail("seek failed");
    long size=ftell(f);
    if(size<0 || size>128*1024*1024)fail("invalid experiment input length");
    *n=(size_t)size;rewind(f);
    unsigned char *p=malloc(*n?*n:1);if(!p || fread(p,1,*n,f)!=*n)fail("read failed");fclose(f);return p;
}
int main(int argc,char **argv) {
    if(argc!=5)fail("usage: pathbench-openzl MODE LEVEL INPUT OUTPUT");
    mode=atoi(argv[1]);level=atoi(argv[2]);if(mode<0 || mode>9)fail("invalid mode");
    size_t n;unsigned char *input=read_file(argv[3],&n);
    empty=n==0;
    if(mode!=3 && n%8)fail("typed input must contain little-endian u64 values");
    size_t capacity=mode==3?ZSTD_compressBound(n):ZL_compressBound(n);
    unsigned char *out=malloc(capacity?capacity:1),*decoded=malloc(n?n:1);if(!out || !decoded)fail("allocation failed");
    size_t size=0;double start=now();
    for(int i=0;i<3;i++) {
        if(mode==3) {size=ZSTD_compress(out,capacity,input,n,level);if(ZSTD_isError(size))fail("zstd error");}
        else size=checked(ZL_compress_usingGraphFn(out,capacity,input,n,graph));
    }
    double encode_seconds=(now()-start)/3;
    FILE *f=fopen(argv[4],"wbx");if(!f)fail("output exists or cannot be created");
    if(fwrite(out,1,size,f)!=size || fclose(f))fail("write failed");
    free(out);
    size_t persisted;out=read_file(argv[4],&persisted);if(persisted!=size)fail("persisted size mismatch");
    start=now();
    for(int i=0;i<101;i++) {
        size_t got=mode==3?ZSTD_decompress(decoded,n,out,size):checked(ZL_decompress(decoded,n,out,size));
        if(got!=n)fail("decode size mismatch");
    }
    double decode_seconds=(now()-start)/101;
    if(memcmp(decoded,input,n))fail("roundtrip mismatch");
    printf("{\"bytes\":%zu,\"raw_bytes\":%zu,\"encode_seconds\":%.9f,\"decode_seconds\":%.9f,\"verified\":true}\n",size,n,encode_seconds,decode_seconds);
    free(input);free(out);free(decoded);return 0;
}
