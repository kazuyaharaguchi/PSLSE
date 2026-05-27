/*
 * cpu_time() は、プロセスの CPU 消費時間を秒単位で返します。
 *
 * プログラム開始直後に、ダミーの呼び出しを入れる方が安全です。
 * 関数の型は double なので、小数付きの値が返ります。
 * エラーがあれば、0.0 を返しますが、そんなことはまず起こらないでしょう。
 * gcc/egcs (Sun,Linux,FreeBSD), Turbo C/LSI-C (MS-DOS) で作動確認済み。
 * 最新版は以下で公開します。
 * http://www.nn.iij4u.or.jp/~tutimura/c/
 * 2006/ 1/14 (C) tutimura(a)nn.iij4u.or.jp (^^) ご自由にお使い下さい。
 */

/* 分割コンパイルしない時の使い方
 *
 *  #include "cpu_time.c"
 *  ...
 *  main() {
 *      cpu_time(); <== これはダミーの呼出し
 *      ...
 *      printf( "time:%.2f(sec)\n", cpu_time() );
 *  }
 */

/* 分割コンパイルする時の使い方
 *
 * 共通ヘッダに以下を記述する。
 *  double cpu_time( void );
 */


#ifdef TEST
#include <stdio.h>
#define TEST_MESSAGE(s) printf("by %s\n",s)
#else
#define TEST_MESSAGE(s)
#endif /* TEST */

#ifdef __STDC__
double cpu_time( void );
#endif /* __STDC__ */


#if ( (defined(MSDOS) || defined(__MSDOS__) || defined(__BORLANDC__) || \
       defined(LSI_C) || defined(_MSC_VER)  || defined(__WATCOMC__)) && \
      !(defined(unix) || defined(__unix))  )

	/* ↓MS-DOS↓ */
#include <time.h>
#if !defined(CLOCKS_PER_SEC) && defined(CLK_TCK)
#define CLOCKS_PER_SEC CLK_TCK
#endif
double cpu_time( void ) {
    clock_t t;  static clock_t last = (clock_t)-1;
    TEST_MESSAGE("clock()");
    t = clock();
    if (last == (clock_t)-1) last = t;
    return (double)(t-last)/CLOCKS_PER_SEC;
}
	/* ↑MS-DOS↑ */

#else /* unix */

#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h>
#ifdef RUSAGE_SELF

	/* ↓Sun 5.5, Linux, FreeBSD, DJGPP↓ */
double cpu_time() {
    struct rusage tmp;
    TEST_MESSAGE("getrusage()");
    if ( getrusage( RUSAGE_SELF, &tmp ) ) return 0.0;
    return (double)(tmp.ru_utime.tv_sec)
          +(double)(tmp.ru_utime.tv_usec)/1000000;
}
	/* ↑Sun 5.5, Linux, FreeBSD, DJGPP↑ */
#else
	/* ↓Sun 5.3↓ */
#include <sys/types.h>
#include <sys/times.h>
double cpu_time() {
    struct tms tmp;
    TEST_MESSAGE("times()");
    if ( times( &tmp ) == -1 ) return 0.0;
    return (double)(tmp.tms_utime)/CLK_TCK;
}
	/* ↑Sun 5.3↑ */
#endif

#endif /* unix */


#ifdef TEST
int main() {
    int i, j, k;

    cpu_time(); /* dummy */
    for ( i=0; i<25; i++ ) {
        for ( j=1; j<300000; j++ ) {
            k = j*j*j*j;
            k /= j;
            k /= j;
            k /= j;
        }
        j=k;
        printf( "time:%.4f(sec)\n", cpu_time() );
    }
    return 0;
}
#endif /* TEST */


/* 詳細説明
 *
 * このプログラムは、ANSI 規格で
 *     (double) clock()/CLOCKS_PER_SEC
 * に相当する計算を行いますが、
 * 以下のような不具合を解消しています。
 *
 *  -- 古いANSI規格では CLOCKS_PER_SEC ではなく CLK_TCK が
 *     定義されていた。
 *  --（主にUnix系の）処理系によっては CLOCKS_PER_SEC が
 *     1000000 と定義されていて、 long(32bit)の clock() が
 *     30分ちょっとで桁溢れを起こす。
 *  -- CPU timeではなく、現在時刻を返す処理系がある
 *     （プログラム開始直後に0でない）。
 *     半面、初めての呼び出しでは（開始直後でなくても）
 *     必ず 0 を返す処理系もある。
 */

/* 動作確認環境
 *
 * SUN 系
 * gcc 2.5.8 (SunOS 4.1.3) -- getrusage() で警告あり
 * gcc 2.6.3 (SunOS 5.5/5.5.1) -- getrusage() で警告あり
 * egcs 1.1.2 (SunOS 5.5/5.5.1)
 *
 * Linux 系
 * gcc 2.7.2.3 (Linux 2.0/2.2, glibc2.0.7)
 * egcs 1.0.3 (Linux 2.0/2.2, glibc2.0.7)
 *
 * FreeBSD 系
 * gcc 2.7.2.1 (FreeBSD 3.2-RELEASE)
 *
 * MS-DOS 系
 * DJGPP
 * Turbo C 2nd Ver.1.0
 * Turbo C++ Ver.2.0
 * Turbo C++ Ver.4.0J
 * LSI-C86 Ver.3.3c
 *
 * Windows 系
 * Visual C++ Ver 6.0
 * Borland C compiler 5.5
 */

/* 歴史
 * 1997/ 2/20  my_time() として登場。
 * 2000/ 2/ 6  cpu_time() に改名。
 *             大幅なスタイル変更するも、処理内容はそのまま。
 * 2002/ 3/11  MS-DOS 部分で初回呼び出しを見分けるために
 *             last に 0 を代入していたが、よく考えると
 *             これでは初回と見分けがつかないことがあるので
 *             -1 を代入することに変更した。
 * 2002/11/29  Visual C++ に対応した。
 * 2006/ 1/10  Borland C++, WATCOM C/C++ 用のマクロを書き加えた。
 *             テストルーチンのループ回数を増やして、割算も追加した。
 * 2006/ 1/14  Borland C++ の動作報告をいただいた。
 */
