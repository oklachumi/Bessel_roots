!*********************************************************************
!*   Compute the zeros of Bessel functions Jn(x), Yn(x), and their   *
!*   derivatives using subroutine JYZO                               *
!* ----------------------------------------------------------------- *
!* SAMPLE RUN:                                                       *
!* (Compute 10 zeroes for n=1).                                      *
!*                                                                   *
!*  Please enter order n and number of zeroes: 1 10                  *
!*                                                                   *
!*  Zeros of Bessel functions Jn(x), Yn(x) and their derivatives     *
!*                       ( n = 1 )                                   *
!*   m       jnm           j'nm          ynm           y'nm          *
!*  -----------------------------------------------------------      *
!*   1     3.8317060     1.8411838     2.1971413     3.6830229       *
!*   2     7.0155867     5.3314428     5.4296810     6.9415000       *
!*   3    10.1734681     8.5363164     8.5960059    10.1234047       *
!*   4    13.3236919    11.7060049    11.7491548    13.2857582       *
!*   5    16.4706301    14.8635886    14.8974421    16.4400580       *
!*   6    19.6158585    18.0155279    18.0434023    19.5902418       *
!*   7    22.7600844    21.1643699    21.1880689    22.7380347       *
!*   8    25.9036721    24.3113269    24.3319426    25.8843146       *
!*   9    29.0468285    27.4570506    27.4752950    29.0295758       *
!*  10    32.1896799    30.6019230    30.6182865    32.1741182       *
!*  -----------------------------------------------------------      *
!*                                                                   *
!* ----------------------------------------------------------------- *
!* Ref.: www.esrf.fr/computing/expg/libraries/smf/PROGRAMS/MJYZO.FOR *
!*                                                                   *
!*                             F90 Release 1.0 By J-P Moreau, Paris. *
!*                                      (www.jpmoreau.fr)            *
!*********************************************************************
program mjyzo
implicit real*8(a-h,j-z)
dimension rj0(101),rj1(101),ry0(101),ry1(101)

open(unit=78,file='result.txt')
write(*,10,advance ='no')
10  format('please enter order n and nth zero roots: ')
read(*,*)in,int
write(*,*)

write(*,'(9(1x,a14))')'n','m','jnm','jnm_deri','ynm','ynm_deri'   ! deri = derivative
write(78,'(9(1x,a14))')'n','m','jnm','jnm_deri','ynm','ynm_deri'
do iin=0,in
  call jyzo(iin,int,rj0,rj1,ry0,ry1)
  do im=1,int
    write(*,'(2(1x,i14),4(1x,1pd14.7))')iin,im,rj0(im),rj1(im),ry0(im),ry1(im)
    write(78,'(2(1x,i14),4(1x,1pd14.7))')iin,im,rj0(im),rj1(im),ry0(im),ry1(im)
  enddo
enddo

stop
end

subroutine jyzo(n,nt,rj0,rj1,ry0,ry1)
!       ======================================================

!       Purpose: Compute the zeros of Bessel functions Jn(x),

!                Yn(x), and their derivatives

!       Input :  n  --- Order of Bessel functions (0 to 100)

!                NT --- Number of zeros (roots)

!       Output:  RJ0(L) --- L-th zero of Jn(x),  L=1,2,...,NT

!                RJ1(L) --- L-th zero of Jn'(x), L=1,2,...,NT

!                RY0(L) --- L-th zero of Yn(x),  L=1,2,...,NT

!                RY1(L) --- L-th zero of Yn'(x), L=1,2,...,NT

!       Routine called: JYNDD for computing Jn(x), Yn(x), and

!                       their first and second derivatives

!       ======================================================

implicit double precision (a-h,o-z)
dimension rj0(nt),rj1(nt),ry0(nt),ry1(nt)

if (n.le.20) then
    x=2.82141+1.15859*n
else
    x=n+1.85576*n**0.33333+1.03315/n**0.33333
endif

l=0
10 x0=x
call jyndd(n,x,bjn,djn,fjn,byn,dyn,fyn)
x=x-bjn/djn
if (dabs(x-x0).gt.1.0d-9) go to 10
l=l+1
rj0(l)=x
x=x+3.1416+(0.0972+0.0679*n-0.000354*n**2)/l
if (l.lt.nt) go to 10

if (n.le.20) then
    x=0.961587+1.07703*n
else
    x=n+0.80861*n**0.33333+0.07249/n**0.33333
endif

if (n.eq.0) x=3.8317
l=0
15 x0=x
call jyndd(n,x,bjn,djn,fjn,byn,dyn,fyn)
x=x-djn/fjn
if (dabs(x-x0).gt.1.0d-9) go to 15
l=l+1
rj1(l)=x
x=x+3.1416+(0.4955+0.0915*n-0.000435*n**2)/l
if (l.lt.nt) go to 15

if (n.le.20) then
    x=1.19477+1.08933*n
else
    x=n+0.93158*n**0.33333+0.26035/n**0.33333
endif
           
l=0
20 x0=x
call jyndd(n,x,bjn,djn,fjn,byn,dyn,fyn)
x=x-byn/dyn
if (dabs(x-x0).gt.1.0d-9) go to 20
l=l+1
ry0(l)=x
x=x+3.1416+(0.312+0.0852*n-0.000403*n**2)/l
if (l.lt.nt) go to 20

if (n.le.20) then
    x=2.67257+1.16099*n
else
    x=n+1.8211*n**0.33333+0.94001/n**0.33333
endif
 
l=0
25 x0=x
call jyndd(n,x,bjn,djn,fjn,byn,dyn,fyn)
x=x-dyn/fyn
if (dabs(x-x0).gt.1.0d-9) go to 25
l=l+1
ry1(l)=x
x=x+3.1416+(0.197+0.0643*n-0.000286*n**2)/l 
if (l.lt.nt) go to 25

return
end

subroutine jyndd(n,x,bjn,djn,fjn,byn,dyn,fyn)
!       ===========================================================

!       Purpose: Compute Bessel functions Jn(x) and Yn(x), and

!                their first and second derivatives 

!       Input:   x   ---  Argument of Jn(x) and Yn(x) ( x > 0 )

!                n   ---  Order of Jn(x) and Yn(x)

!       Output:  BJN ---  Jn(x)

!                DJN ---  Jn'(x)

!                FJN ---  Jn"(x)

!                BYN ---  Yn(x)

!                DYN ---  Yn'(x)

!                FYN ---  Yn"(x)

!       ===========================================================

implicit double precision (a-h,o-z)
dimension bj(102),by(102)

do nt=1,900
    mt=int(0.5*log10(6.28*nt)-nt*log10(1.36*dabs(x)/nt))
    if (mt.gt.20) go to 15
enddo
15 m=nt
bs=0.0d0
f0=0.0d0
f1=1.0d-35
su=0.0d0

do 20 k=m,0,-1
    f=2.0d0*(k+1.0d0)*f1/x-f0
    if (k.le.n+1) bj(k+1)=f
    if (k.eq.2*int(k/2)) then
        bs=bs+2.0d0*f
        if (k.ne.0) su=su+(-1)**(k/2)*f/k
    endif
    f0=f1
20 f1=f

do 25 k=0,n+1
25 bj(k+1)=bj(k+1)/(bs-f)
bjn=bj(n+1)
ec=0.5772156649015329d0
e0=0.3183098861837907d0
s1=2.0d0*e0*(dlog(x/2.0d0)+ec)*bj(1)
f0=s1-8.0d0*e0*su/(bs-f)
f1=(bj(2)*f0-2.0d0*e0/x)/bj(1)
by(1)=f0
by(2)=f1

do 30 k=2,n+1
    f=2.0d0*(k-1.0d0)*f1/x-f0
    by(k+1)=f
    f0=f1
30 f1=f
byn=by(n+1)
djn=-bj(n+2)+n*bj(n+1)/x
dyn=-by(n+2)+n*by(n+1)/x
fjn=(n*n/(x*x)-1.0d0)*bjn-djn/x
fyn=(n*n/(x*x)-1.0d0)*byn-dyn/x

return
end
