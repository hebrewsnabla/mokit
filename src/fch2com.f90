! written by jxzou at 20201210: copy from fch2inporb.f90 and do modifications
!  generate Molpro .com file from Gaussian .fch(k) file
! updated by jxzou at 20210407: remove '-uhf' and add automatic determination

program main
 use fch_content, only: iout
 use util_wrapper, only: fch2inp_wrap
 implicit none
 integer :: i, system
 character(len=240) :: fchname, inpname
 logical :: sph

 i = iargc()
 if(i /= 1) then
  write(iout,'(/,A)') ' ERROR in subroutine fch2com: wrong command line argument!'
  write(iout,'(A,/)') ' Example (R(O)HF, GVB, CAS): fch2com a.fch'
  stop
 end if

 fchname = ' '
 sph = .true.

 call getarg(1, fchname)
 i = index(fchname, '.fch', back=.true.)
 if(i == 0) then
  write(iout,'(A)') "ERROR in subroutine fch2com: no '.fch' suffix in&
                   & filename="//TRIM(fchname)
  stop
 end if
 inpname = fchname(1:i-1)//'.inp'

 call require_file_exist(fchname)
 call fch2inp_wrap(fchname, .false., 0, 0) ! generate GAMESS .inp file

 call check_sph(fchname, sph)
 if(sph) then
  i = system('bas_gms2molpro '//TRIM(inpname)//' -sph')
 else
  i = system('bas_gms2molpro '//TRIM(inpname))
 end if

 if(i /= 0) then
  write(iout,'(/,A)') 'ERROR in subroutine fch2com: failed to call utility&
                     & bas_gms2molpro. Two possible reasons:'
  write(iout,'(A)')   '(1) The file '//TRIM(fchname)//' may be incomplete.'
  write(iout,'(A,/)') '(2) You forgot to compile the utility bas_gms2molpro.'
  stop
 end if

 call delete_file(inpname)
 call fch2com(fchname)
 stop
end program main

! nbf: the number of basis functions
! nif: the number of independent functions, i.e., the number of MOs

! read the MOs in .fch(k) file and adjust its d,f,g, etc. functions order
!  of Gaussian to that of Molcas
subroutine fch2com(fchname)
 use fch_content, only: iout, check_uhf_in_fch
 implicit none
 integer :: i, j, k, length, orbid
 integer :: nalpha, nbeta, nbf, nif
 integer :: nbf0, nbf1
 integer :: n10fmark, n15gmark, n21hmark
 integer :: n5dmark, n7fmark, n9gmark, n11hmark
 integer, allocatable :: shell_type(:), shell2atom_map(:)
 ! mark the index where d, f, g, h functions begin
 integer, allocatable :: d_mark(:), f_mark(:), g_mark(:), h_mark(:)
 character(len=240) :: fileA, fileB
 character(len=240), intent(in) :: fchname
 real(kind=8), allocatable :: coeff(:,:)
 logical :: uhf, sph

 call read_na_and_nb_from_fch(fchname, nalpha, nbeta)
 call read_nbf_and_nif_from_fch(fchname, nbf, nif)
 nbf0 = nbf

 ! read MO Coefficients
 call check_uhf_in_fch(fchname, uhf)
 if(uhf) then
  allocate(coeff(nbf,2*nif))
  call read_mo_from_fch(fchname, nbf, nif, 'a', coeff(:,1:nif))
  call read_mo_from_fch(fchname, nbf, nif, 'b', coeff(:,nif+1:2*nif))
  nif = 2*nif   ! double the size
 else
  allocate(coeff(nbf,nif))
  call read_mo_from_fch(fchname, nbf, nif, 'a', coeff)
 end if

 call read_ncontr_from_fch(fchname, k)
 allocate(shell_type(2*k), source=0)
 allocate(shell2atom_map(2*k), source=0)
 call read_shltyp_and_shl2atm_from_fch(fchname, k, shell_type, shell2atom_map)

 if(ANY(shell_type<-1) .and. ANY(shell_type>1)) then
  write(iout,'(A)') 'ERROR in subroutine fch2com: mixed spherical harmonic/&
                   &Cartesian functions detected.'
  write(iout,'(A)') 'You probably used a basis set like 6-31G(d) in Gaussian. Its&
                   & default setting is (6D,7F).'
  write(iout,'(A)') "You need to add '5D 7F' or '6D 10F' keywords in Gaussian."
  stop
 else if( ANY(shell_type<-1) ) then
  sph = .true.
 else
  sph = .false.
 end if

! first we adjust the basis functions in each MO according to the Shell to atom map
 ! 1) split the 'L' into 'S' and 'P', this is to ensure that D comes after L functions
 call split_L_func(k, shell_type, shell2atom_map, length)

 ! 2) sort the shell_type and shell2atom_map by ascending order
 ! MOs will be adjusted accordingly
 call sort_shell_and_mo(length, shell_type, shell2atom_map, nbf, nif, coeff)
! adjust done

 deallocate(shell2atom_map)

! then we adjust the basis functions in each MO according to the type of basis functions
 k = length  ! update k
 n10fmark = 0
 n15gmark = 0
 n21hmark = 0
 n5dmark = 0
 n7fmark = 0
 n9gmark = 0
 n11hmark = 0
 allocate(d_mark(k), f_mark(k), g_mark(k), h_mark(k))
 d_mark = 0
 f_mark = 0
 g_mark = 0
 h_mark = 0
 nbf = 0
 do i = 1, k, 1
  select case(shell_type(i))
  case( 0)   ! S
   nbf = nbf + 1
  case( 1)   ! 3P
   nbf = nbf + 3
  case(-1)   ! SP or L
   nbf = nbf + 4
  case(-2)   ! 5D
   n5dmark = n5dmark + 1
   d_mark(n5dmark) = nbf + 1
   nbf = nbf + 5
  case( 2)   ! 6D
   nbf = nbf + 6
  case(-3)   ! 7F
   n7fmark = n7fmark + 1
   f_mark(n7fmark) = nbf + 1
   nbf = nbf + 7
  case( 3)   ! 10F
   n10fmark = n10fmark + 1
   f_mark(n10fmark) = nbf + 1
   nbf = nbf + 10
  case(-4)   ! 9G
   n9gmark = n9gmark + 1
   g_mark(n9gmark) = nbf + 1
   nbf = nbf + 9
  case( 4)   ! 15G
   n15gmark = n15gmark + 1
   g_mark(n15gmark) = nbf + 1
   nbf = nbf + 15
  case(-5)   ! 11H
   n11hmark = n11hmark + 1
   h_mark(n11hmark) = nbf + 1
   nbf = nbf + 11
  case( 5)   ! 21H
   n21hmark = n21hmark + 1
   h_mark(n21hmark) = nbf + 1
   nbf = nbf + 21
  case default
   write(iout,'(A)') 'ERROR in subroutine fch2com: shell_type(i) out of range.'
   write(iout,'(2(A,I0))') 'i=', i, ', k=', k
   stop
  end select
 end do ! for i
 deallocate(shell_type)

 ! adjust the order of d, f, etc. functions
 do i = 1, n5dmark, 1
  call fch2com_permute_5d(nif,coeff(d_mark(i):d_mark(i)+4,:))
 end do
 do i = 1, n7fmark, 1
  call fch2com_permute_7f(nif,coeff(f_mark(i):f_mark(i)+6,:))
 end do
 do i = 1, n10fmark, 1
  call fch2com_permute_10f(nif,coeff(f_mark(i):f_mark(i)+9,:))
 end do
 do i = 1, n9gmark, 1
  call fch2com_permute_9g(nif,coeff(g_mark(i):g_mark(i)+8,:))
 end do
 do i = 1, n15gmark, 1
  call fch2com_permute_15g(nif,coeff(g_mark(i):g_mark(i)+14,:))
 end do
 do i = 1, n11hmark, 1
  call fch2com_permute_11h(nif,coeff(h_mark(i):h_mark(i)+10,:))
 end do
 do i = 1, n21hmark, 1
  call fch2com_permute_21h(nif,coeff(h_mark(i):h_mark(i)+20,:))
 end do
! adjustment finished

 deallocate(d_mark, f_mark, g_mark, h_mark)

! print MOs into a plain text file
 fileA = fchname
 fileB = fchname
 call convert2molpro_fname(fileA, '.a')
 call convert2molpro_fname(fileB, '.b')

 if(uhf) nif = nif/2

 open(newunit=orbid,file=TRIM(fileA),status='replace')
 do i = 1, nbf0, 1
  write(orbid,'(5(1X,ES16.9))') (coeff(i,j),j=1,nif)
 end do ! for i
 close(orbid)

 if(uhf) then
  open(newunit=orbid,file=TRIM(fileB),status='replace')
  do i = 1, nbf0, 1
   write(orbid,'(5(1X,ES16.9))') (coeff(i,j),j=nif+1,2*nif)
  end do ! for i
  close(orbid)
 end if
! print done

 deallocate(coeff)
 return
end subroutine fch2com

! split the 'L' into 'S' and 'P'
subroutine split_L_func(k, shell_type, shell2atom_map, length)
 implicit none
 integer i, k0
 integer,intent(in) :: k
 integer,intent(inout) :: shell_type(2*k), shell2atom_map(2*k)
 integer,intent(out) :: length
 integer,allocatable :: temp1(:), temp2(:)

 k0 = 2*k
 length = k
 ! set initial values for arrays shell_type, assume 15 has not be used
 shell_type(k+1:k0) = 15
 i = 1
 do while(shell_type(i) /= 15)
  if(shell_type(i) /= -1) then
   i = i + 1
   cycle
  end if
  shell_type(i) = 0
  allocate( temp1(i+1 : k0-1), temp2(i+1 : k0-1) )
  temp1(i+1 : k0-1) = shell_type(i+1 : k0-1)
  shell_type(i+2 : k0) = temp1(i+1 : k0-1)
  temp2(i+1 : k0-1) = shell2atom_map(i+1 : k0-1)
  shell2atom_map(i+2 : k0) = temp2(i+1 : k0-1)
  deallocate(temp1, temp2)
  shell_type(i+1) = 1
  shell2atom_map(i+1) = shell2atom_map(i)
  i = i + 2
 end do

 length = i - 1
 shell_type(i : k0) = 0
 return
end subroutine split_L_func

! sort the shell_type, shell2atom_map by ascending order
! MOs will be adjusted accordingly
subroutine sort_shell_and_mo(ilen, shell_type, shell2atom_map, nbf, nif, coeff2)
 implicit none
 integer :: i, j, k
 integer :: ibegin, iend, natom
 integer :: jbegin, jend
 integer, parameter :: ntype = 10
 integer, parameter :: num0(ntype) = [0, 1, -2, 2, -3, 3, -4, 4, -5, 5]
 integer, parameter :: num1(ntype) = [1, 3, 5, 6, 7, 10, 9, 15, 11, 21]
 !                                    S  P  5D 6D 7F 10F 9G 15G 11H 21H
 integer num(ntype)

 integer, intent(in) :: ilen, nbf, nif
 integer, intent(inout) :: shell_type(ilen), shell2atom_map(ilen)
 integer, allocatable :: ith(:), ith_bas(:), tmp_type(:)
 real(kind=8), intent(inout) :: coeff2(nbf,nif)

 ! find the number of atoms
 natom = shell2atom_map(ilen)

 allocate(ith(0:natom), ith_bas(0:natom))
 ith = 0
 ith_bas = 0

 ! find the end position of each atom in array shell2atom_map
 do i = 1, natom, 1
  ith(i) = count(shell2atom_map==i) + ith(i-1)
 end do

 ! find the end position of basis functions between two atoms
 do i = 1, natom, 1
  ibegin = ith(i-1) + 1
  iend = ith(i)
  allocate(tmp_type(ibegin:iend))
  tmp_type = 0
  tmp_type = shell_type(ibegin:iend)
  num = 0
  do j = 1, ntype, 1
   num(j) = count(tmp_type == num0(j))
  end do
  k = 0
  do j = 1, ntype, 1
   k = k + num(j)*num1(j)
  end do
  ith_bas(i) = ith_bas(i-1) + k
  deallocate(tmp_type)
 end do

 ! adjust the MOs in each atom
 do i = 1, natom, 1
  ibegin = ith(i-1) + 1
  iend = ith(i)
  jbegin = ith_bas(i-1) + 1
  jend = ith_bas(i)
  call sort_shell_and_mo_in_each_atom(iend-ibegin+1, shell_type(ibegin:iend), &
  & jend-jbegin+1, nif, coeff2(jbegin:jend,1:nif))
 end do
 ! adjust the MOs in each atom done
 deallocate(ith, ith_bas)
 return
end subroutine sort_shell_and_mo

subroutine sort_shell_and_mo_in_each_atom(ilen1, shell_type, ilen2, nif, coeff2)
 implicit none
 integer i, tmp_type
 integer ibegin, iend, jbegin, jend
 integer,parameter :: ntype = 10
 integer,parameter :: num0(ntype) = [0, 1, -2, 2, -3, 3, -4, 4, -5, 5]
 integer,parameter :: num1(ntype) = [1, 3, 5, 6, 7, 10, 9, 15, 11, 21]
 !                                   S  P  5D 6D 7F 10F 9G 15G 11H 21H
 integer,parameter :: rnum(-5:5) = [9, 7, 5, 3, 0, 1, 2, 4, 6, 8, 10]

 integer,intent(in) :: nif, ilen1, ilen2
 integer,intent(inout) :: shell_type(ilen1)
 integer,allocatable :: ith_bas(:)
 real(kind=8),intent(inout) :: coeff2(ilen2,nif)
 real(kind=8),allocatable :: tmp_coeff1(:,:), tmp_coeff2(:,:)
 logical sort_done

 tmp_type = 0

 ! find the end position of basis functions within an atom
 allocate(ith_bas(0:ilen1))
 ith_bas = 0
 do i = 1, ilen1, 1
   ith_bas(i) = ith_bas(i-1) + num1(rnum(shell_type(i)))
 end do

 sort_done = .false.
 do while(.not. sort_done)
  sort_done = .true.
  do i = 1, ilen1-1, 1
    if(shell_type(i) == 0) cycle
    if(ABS(shell_type(i+1)) >= ABS(shell_type(i))) cycle
    sort_done = .false.
    tmp_type = shell_type(i+1)
    shell_type(i+1) = shell_type(i)
    shell_type(i) = tmp_type
    ibegin = ith_bas(i-1) + 1
    iend = ith_bas(i)
    jbegin = ith_bas(i) + 1
    jend = ith_bas(i+1)
    allocate(tmp_coeff1(ibegin:iend,nif), tmp_coeff2(jbegin:jend,nif))
    tmp_coeff1 = 0.0d0
    tmp_coeff2 = 0.0d0
    tmp_coeff1 = coeff2(ibegin:iend,:)
    tmp_coeff2 = coeff2(jbegin:jend,:)
    ith_bas(i) = ibegin+jend-jbegin
    coeff2(ibegin: ith_bas(i),:) = tmp_coeff2
    ith_bas(i+1) = jend+iend-jbegin+1
    coeff2(ibegin+jend-jbegin+1: ith_bas(i+1),:) = tmp_coeff1
    deallocate(tmp_coeff1, tmp_coeff2)
  end do
 end do
 deallocate(ith_bas)
 return
end subroutine sort_shell_and_mo_in_each_atom

subroutine fch2com_permute_5d(nif,coeff)
 implicit none
 integer :: i
 integer, parameter :: order(5) = [1,5,2,4,3]
 integer, intent(in) :: nif
 real(kind=8), intent(inout) :: coeff(5,nif)
 real(kind=8), allocatable :: coeff2(:,:)
! From: the order of spherical d functions in Gaussian
! To: the order of spherical d functions in Molpro
! 1    2    3    4    5
! d0 , d+1, d-1, d+2, d-2
! d0 , d2-, d1+, d2+, d1-

 allocate(coeff2(5,nif), source=0.0d0)
 forall(i = 1:5) coeff2(i,:) = coeff(order(i),:)
 coeff = coeff2
 deallocate(coeff2)
 return
end subroutine fch2com_permute_5d

subroutine fch2com_permute_7f(nif,coeff)
 implicit none
 integer :: i
 integer, parameter :: order(7) = [2,3,1,6,5,7,4]
 integer, intent(in) :: nif
 real(kind=8), intent(inout) :: coeff(7,nif)
 real(kind=8), allocatable :: coeff2(:,:)
! From: the order of spherical f functions in Gaussian
! To: the order of spherical f functions in Molpro
! 1    2    3    4    5    6    7
! f0 , f+1, f-1, f+2, f-2, f+3, f-3
! f1+, f1-, f0 , f3+, f2-, f3-, f2+

 allocate(coeff2(7,nif), source=0.0d0)
 forall(i = 1:7) coeff2(i,:) = coeff(order(i),:)
 coeff = coeff2
 deallocate(coeff2)
 return
end subroutine fch2com_permute_7f

subroutine fch2com_permute_10f(nif,coeff)
 implicit none
 integer :: i
 integer, parameter :: order(10) = [1,2,3,5,6,4,9,7,8,10]
 integer, intent(in) :: nif
 real(kind=8), intent(inout) :: coeff(10,nif)
 real(kind=8), allocatable :: coeff2(:,:)
! From: the order of Cartesian f functions in Gaussian
! To: the order of Cartesian f functions in Molpro
! 1   2   3   4   5   6   7   8   9   10
! XXX,YYY,ZZZ,XYY,XXY,XXZ,XZZ,YZZ,YYZ,XYZ
! xxx,yyy,zzz,xxy,xxz,xyy,yyz,xzz,yzz,xyz

 allocate(coeff2(10,nif), source=0.0d0)
 forall(i = 1:10) coeff2(i,:) = coeff(order(i),:)
 coeff = coeff2
 deallocate(coeff2)
 return
end subroutine fch2com_permute_10f

subroutine fch2com_permute_9g(nif,coeff)
 implicit none
 integer :: i
 integer, parameter :: order(9) = [1,5,2,8,3,4,9,6,7]
 integer, intent(in) :: nif
 real(kind=8), intent(inout) :: coeff(9,nif)
 real(kind=8), allocatable :: coeff2(:,:)
! From: the order of spherical g functions in Gaussian
! To: the order of spherical g functions in Molpro
! 1    2    3    4    5    6    7    8    9
! g0 , g+1, g-1, g+2, g-2, g+3, g-3, g+4, g-4
! g0 , g2-, g1+, g4+, g1-, g2+, g4-, g3+, g3-

 allocate(coeff2(9,nif), source=0.0d0)
 forall(i = 1:9) coeff2(i,:) = coeff(order(i),:)
 coeff = coeff2
 deallocate(coeff2)
 return
end subroutine fch2com_permute_9g

subroutine fch2com_permute_15g(nif,coeff)
 implicit none
 integer :: i
 integer, parameter :: order(15) = [15,5,1,14,13,9,4,6,2,12,10,3,11,8,7]
 integer, intent(in) :: nif
 real(kind=8), intent(inout) :: coeff(15,nif)
 real(kind=8), allocatable :: coeff2(:,:)
! From: the order of Cartesian g functions in Gaussian
! To: the order of Cartesian g functions in Molpro
! 1    2    3    4    5    6    7    8    9    10   11   12   13   14   15
! ZZZZ,YZZZ,YYZZ,YYYZ,YYYY,XZZZ,XYZZ,XYYZ,XYYY,XXZZ,XXYZ,XXYY,XXXZ,XXXY,XXXX
! xxxx,yyyy,zzzz,xxxy,xxxz,xyyy,yyyz,xzzz,yzzz,xxyy,xxzz,yyzz,xxyz,xyyz,xyzz

 allocate(coeff2(15,nif), source=0.0d0)
 forall(i = 1:15) coeff2(i,:) = coeff(order(i),:)
 coeff = coeff2
 deallocate(coeff2)
 return
end subroutine fch2com_permute_15g

subroutine fch2com_permute_11h(nif,coeff)
 implicit none
 integer :: i
 integer, parameter :: order(11) = [2,3,4,6,9,7,8,11,1,10,5]
 integer, intent(in) :: nif
 real(kind=8), intent(inout) :: coeff(11,nif)
 real(kind=8), allocatable :: coeff2(:,:)
! From: the order of spherical h functions in Gaussian
! To: the order of spherical h functions in Molpro
! 1    2    3    4    5    6    7    8    9    10   11
! h0 , h+1, h-1, h+2, h-2, h+3, h-3, h+4, h-4, h+5, h-5
! h1+, h1-, h2+, h3+, h4-, h3-, h4+, h5-, h0 , h5+, h2-

 allocate(coeff2(11,nif), source=0.0d0)
 forall(i = 1:11) coeff2(i,:) = coeff(order(i),:)
 coeff = coeff2
 deallocate(coeff2)
 return
end subroutine fch2com_permute_11h

subroutine fch2com_permute_21h(nif,coeff)
 implicit none
 integer :: i
 integer, intent(in) :: nif
 real(kind=8), intent(inout) :: coeff(21,nif)
 real(kind=8), allocatable :: coeff2(:,:)
! From: the order of Cartesian h functions in Gaussian
! To: the order of Cartesian h functions in Molpro
! 1     2     3     4     5     6     7     8     9     10    11    12    13    14    15    16    17    18    19    20    21
! ZZZZZ,YZZZZ,YYZZZ,YYYZZ,YYYYZ,YYYYY,XZZZZ,XYZZZ,XYYZZ,XYYYZ,XYYYY,XXZZZ,XXYZZ,XXYYZ,XXYYY,XXXZZ,XXXYZ,XXXYY,XXXXZ,XXXXY,XXXXX
! xxxxx,xxxxy,xxxxz,xxxyy,xxxyz,xxxzz,xxyyy,xxyyz,xxyzz,xxzzz,xyyyy,xyyyz,xyyzz,xyzzz,xzzzz,yyyyy,yyyyz,yyyzz,yyzzz,yzzzz,zzzzz

 allocate(coeff2(21,nif), source=0.0d0)
 forall(i = 1:21) coeff2(i,:) = coeff(22-i,:)
 coeff = coeff2
 deallocate(coeff2)
 return
end subroutine fch2com_permute_21h

