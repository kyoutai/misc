#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <math.h>
/* 引数解析用 */
typedef struct{
  char infile[1024];
  char outfile[1024];
  int  skip;
}PARAM; 


void Error(){
  printf("Usage: dump2analysis hoge options hoge.lammpstrj\n");
  printf("Options:\n");
  printf(" -h --help        show this help message\n");
  printf("--------------------------------------------------------\n");
  printf(" -s --skip=%%d    Step number\n");
  printf("--------------------------------------------------------\n");
  printf(" -i --in=str      Input file name\n");
  printf(" -o --out=str     Output file name\n");
  printf("--------------------------------------------------------\n");
  exit(1);
}

int SetArgment(int argc, char *argv[], PARAM *arg){
  int opt, option_index;
  
  struct option long_options[] = {
    {"help",  no_argument,  NULL, 'h'},
    {"mode",  required_argument,  NULL, 'm'},
    {"skip",  required_argument,  NULL, 's'},
    {"out",  required_argument,  NULL, 'o'},
    {"in",  required_argument,  NULL, 'i'},
    {0, 0, 0, 0}// 配列の最後はすべて0で埋める
  };
  if (argc == 1) Error(); /* 引数なし */
  /* b,c,d.oは引数が必須 */
  while((opt = getopt_long(argc, argv, "h:s:i:o:", 
                           long_options, &option_index)) != -1){
    switch(opt){
    case 'h': Error(); break;
    case 's': arg->skip = atoi(optarg); break;
    case 'i': strcpy(arg->infile, optarg); break;
    case 'o': strcpy(arg->outfile, optarg); break;
    case '?':
      Error();
    }
  }
  if (arg->skip == 0){
    printf("skip must be > 0!\n");
    Error();
  }
  if (strcmp(arg->outfile, "")==0) {
    printf("outfile (-o) is missing!\n");
    Error();
  }
  if (strcmp(arg->infile, "")==0) {
    printf("infile (-i) is missing!\n");
    Error();
  }
  return 0; 
}
  
int main(int argn, char** argv){
  int i, rstep=0, step=0, atoms;
  char buf[1024];
  FILE *in, *out;
  PARAM param = {"","", 0};
  
  /* 引数の解析 */
  SetArgment(argn, argv, &param);
  /* 原子数の取得 */
  in = fopen(param.infile, "r");
  out = fopen(param.outfile, "w");

  for (i=0; i<3; i++) fgets(buf, sizeof(buf), in);
  fgets(buf, sizeof(buf), in);
  sscanf(buf, "%d", &atoms);
  
  /* 出力 */
  param.skip = param.skip + 1;
  rewind(in);
  while (1){
    if (fgets(buf, sizeof(buf), in) == NULL) break;
    if ((step%param.skip) == 0){
      fprintf(out, "%s", buf);
      for (i=0; i<atoms+8; i++){
    	fgets(buf, sizeof(buf), in);
	fprintf(out, "%s", buf);
      }
      rstep++;
    }
    else {
      for (i=0; i< atoms+8; i++){
    	fgets(buf, sizeof(buf), in);
      }
    }
    step++;
  }
  printf("Total step: %d\n", step);
  printf("Total reduced step: %d\n", rstep);
  
  return 0;
}
