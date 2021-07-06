#include<stdio.h>
#include<stdlib.h>

main()
{
  char a[3];
  FILE *txt;
  int i,j;
  txt=fopen("1.txt","rt");
  for(i=0;;i++) {
  j=fread(a,1,3,txt);
  if(j<3) {printf("j %i\n",j);break;}
  printf("%s %i %i\n",a,i,j);
  }
}
