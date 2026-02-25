#include <stdio.h>

static int prev = -1;
static int first = 1;

void percent(int n,int d,int s, FILE *out)
{
    int x;

    x = (d <= 0 || s <= 0)
      ? 100
      : (int) ( (100.0 * n) / d);
    
    if (n <= 0 || n >= d || x > prev + s){
      prev = x;
      fprintf (out,"%4d%%\b\b\b\b\b",x);
    }

    if (x >= 100){
	fprintf (out,"\n");
	prev = -1;
	first = 1;
    }
    
    return;
}
