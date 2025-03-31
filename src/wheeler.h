
/* file wheeler.h
| copyright (C) 2002-2004 by Robert E. Wheeler
*/

#if !defined(__wheeler_h)
#define __wheeler_h



#define UCHAR unsigned char
#define USHORT unsigned short int
#define UINT unsigned int
#define ULONG unsigned long

#define equals ==
#define repeat do {
#define until(A) }while(!(A))
#define forever }while(1)
#define solongas(A) }while(A)
#define null()

#define true 1
#define false 0


#define maxm(a,b) (((a)>(b))?(a):(b))
#define minm(a,b) (((a)<(b))?(a):(b))
#define absm(a) (((a)<0)?-(a):(a))
#define signm(a) (((a)>0)?1:(((a)<0)?-1:0))
#define SQR(A)   ((A)*(A))

#endif /* Sentinal  */
