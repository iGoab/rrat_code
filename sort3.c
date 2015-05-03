#include <stdio.h>
#include <stdlib.h>

#define M 7
#define NSTACK 50

#define SWAP(a,b) temp=(a);(a)=(b);(b)=temp;

void sort3(int n, float *arr, float *brr, float *crr)
{
	int i, ir, j, jstack, k, l, istack[NSTACK];
	float a, b, c, temp;
	jstack = -1;
	l = 0;
	ir = n - 1;
	for (;;) {
		if ((ir - l) < M) {
			for (j = l + 1; j <= ir; j++) {
				a = arr[j];
				b = brr[j];
				c = crr[j];
				for (i = j - 1; i >=l; i--) {
					if (arr[i] <= a) break;
					arr[i + 1] = arr[i];
					brr[i + 1] = brr[i];
					crr[i + 1] = crr[i];
				}
				i = 0;
				arr[i + 1] = a;
				brr[i + 1] = b;
				crr[i + 1] = c;
			}
			if (jstack < 0) break;
			ir = istack[jstack--];
			l = istack[jstack--];
		} else {
			k = (l + ir) >> 1;
			SWAP(arr[k], arr[l + 1]);
			SWAP(brr[k], brr[l + 1]);
			SWAP(crr[k], crr[l + 1]);
			if (arr[l] > arr[ir]) {
				SWAP(arr[l], arr[ir]);
				SWAP(brr[l], brr[ir]);
				SWAP(crr[l], crr[ir]);
			}
			if (arr[l + 1] > arr[ir]) {
				SWAP(arr[l + 1], arr[ir]);
				SWAP(brr[l + 1], brr[ir]);
				SWAP(crr[l + 1], crr[ir]);
			}
			if (arr[l] > arr[l + 1]) {
				SWAP(arr[l], arr[l + 1]);
				SWAP(brr[l], brr[l + 1]);
				SWAP(crr[l], crr[l + 1]);
			}
			i = l + 1;
			j = ir;
			a = arr[l + 1];
			b = brr[l + 1];
			c = crr[l + 1];
			for (;;) {
				do i++; while (arr[j] < a);
				do j--; while (arr[j] > a);
				if (j < i) break;
				SWAP(arr[i], arr[j]);
				SWAP(brr[i], brr[j]);
				SWAP(crr[i], crr[j]);
			}
			arr[l + 1] = arr[j];
			arr[j] = a;
			brr[l + 1] = brr[j];
			brr[j] = b;
			crr[l + 1] = crr[j];
			crr[j] = c;
			if (jstack >= NSTACK) {
				printf("NSTACK too small in sort3.");
				exit(1);
			}
			if (ir - i + 1 >= j - 1) {
				istack[jstack] = ir;
				istack[jstack - 1] = i;
				ir = j - 1;
			} else {
				istack[jstack] = j - 1;
				istack[jstack - 1] = l;
				l = i;
			}
		}
	} 
}
