program MC

implicit none

integer :: i,i1,i2,j,k,l,pi,pj,n,nb,in,m,ntypes,nid,njumps
integer :: jl,id1,id2,t1,t2,c1,nt1,nt2,totjumps,NP,ii,jj
integer :: ncyc,PF,PT,counter1,hoppf,hoppt,am,hopptn,nvacs
integer :: seed
integer, dimension(8)              ::  rtime
integer, allocatable, dimension(:,:) ::  hist
double precision, allocatable, dimension(:,:) ::  histtot
double precision :: totp,r,E,Beta,T,v0,time,dtime
double precision :: Ea,dE, Einit, Efinal, gamma,Etot,de2ea
!NOTE THAT i IS USED IN A UNCLEVER WAY IN THE MAIN LOOP AND AND i SHOULD NOT BE USED IN ANY SUB-LOOPS!!!

type site
  double precision :: x,y,z
  integer          :: id               ! type of site
  integer          :: species          ! type of species (vacancy , normal oxygen, ad-atom etc.)
  integer          :: nons             ! no of neigbors
  integer          :: njump            ! no of possible jumps
  integer, dimension(200) :: ns ! nieghbohring sites
  integer, dimension(200) :: jt ! jumptype to ns
  double precision, dimension(200)   :: Ea0 ! primary hopping barrier
  double precision, dimension(200)   :: p ! probability of the jump
  double precision, dimension(200)   :: freq ! attempt freq
  integer, dimension(200)            :: jw ! jump where
  integer, dimension(200,4)          :: eff  ! effect of the jump (njump),(site1,site2,newspec1,newspec2)
endtype
 double precision, allocatable,dimension(:,:,:,:,:,:,:) :: barr ! barr(jump_type,id1,id2,type1,type2)
 double precision, allocatable,dimension(:,:,:,:,:,:,:) :: freqs ! barr(jump_type,id1,id2,type1,type2)
 double precision, allocatable,dimension(:,:,:,:,:) :: twobody, amoves ! twobody(jump_type,id1,id2,site1,site2)
 double precision, allocatable,dimension(:,:) :: onebody ! onebody(id,site)
type(site), allocatable, dimension(:) :: sites, isites, fsites


  call DATE_AND_TIME(values=rtime)     ! Get the current time
  seed = rtime(4) * (360000*rtime(5) + 6000*rtime(6) + 100*rtime(7) + rtime(8))
  seed=irand(seed)
  write(*,*), "Random seed:",  seed
  call srand(seed)
!#S# PRINT OUT WELCOME AND INTRUCTIONS
write(*,*)  "                                                         "
write(*,*)  "    - -   - -               - -   - -   - -   - -   - -  "
write(*,*)  "   | | | |           |   | |   | |   | |   | |     |   | "
write(*,*)  "                      - -         - -   - -   - -   - -  "
write(*,*)  "   |   | |           |   | |   | |     |     |     |  \  "
write(*,*)  "          - -               - -               - -        "
write(*,*)  "   Hopping Monte Carlo code by Jolla Kullgren"
write(*,*)  ""
write(*,*)"    There are two main input files: geom.dat and inpu&
t.dat                                                                    "
write(*,*)"                                                     &
                                                                         "
write(*,*)"    -------------------------------------------------&
--------------------------------------------------------------------     "
write(*,*)"     geom.dat:                 Geometry definition   &
                                                                         "
write(*,*)"      Line 1:       i = No of sites                  &
                                                                         "
write(*,*)"      Line 2-i+1:   x y z  type_of_site type_of_speci&
es No_of_neighbors_inc_current_site(N)  N*Id_of_nieghbor N*jump_type     "
write(*,*)"                                                     &
                                                                         "
write(*,*)"      x,x,z are given in any unit e.g. Angstrom.     &
                                                                         "
write(*,*)"      No_of_neighbors is and intereger number giving &
the number of nieghbors including the site itself                        "
write(*,*)"      jump_type  defines how the sites are connected. "
write(*,*)"                                                     &
                                                                         "
write(*,*)"    -------------------------------------------------&
--------------------------------------------------------------------     "
write(*,*)"     input.dat:               2 body energies         &
                                                                         "
write(*,*)"       Line 1:     i = No_of_2_body_energies &
   j=No_of_1_body_energies_moves k=No_of_allow_moves No_of_steps &
    Print_frequency print_type   &
 Temperature  "
write(*,*)"       Line 2:     No_of_different_species  No_of_dif&
ferent_sites No_of_jump_types                                            "
write(*,*)"       Line 3-i+3: jump_type 1st_site 2nd_site  spec_at_1&
 spec_at_2 E "
write(*,*) "      Next j lines:  site species energy"
write(*,*) "      Last k lines:"
write(*,*) "      jump_type 1st_site 2nd_site  spec_at_1&
 spec_at_2 new_spec_at_1 new_spec_at_2 v Ea"
write(*,*)"           1st and 2nd site:  the type of sites invol&
ved in the interaction or allowed moved e.g. surface to bulk            "
write(*,*)"           E: two-body  interaction                           "
write(*,*)"    -------------------------------------------------&
----------------------------------------------------------------------   "

!#E#

!#S# OPEN OUTPUT FILES
open(unit=2,file='KMC.xyz',status='new')
open(unit=3,file='geom.out',status='new')

!#E# OPEN OUTPUT FILES

!#S# READ IN POS
open(unit=10,file='geom.dat',status='old')
read(10,*),n
write(*,*),"    Reading points ",n
allocate(sites(n))
allocate(isites(n))
allocate(fsites(n))


do i=1,n
read(10,*),sites(i)%x, &
sites(i)%y, &
sites(i)%z, &
sites(i)%id, &
sites(i)%species, &
sites(i)%nons, &
(sites(i)%ns(j), j=1,sites(i)%nons), &
(sites(i)%jt(j), j=1,sites(i)%nons)

enddo
close(10)
write(*,*), "    Complete"
write(*,*), "    ---------------------------------------------"


! #E# READ IN POS

! #S# READ IN PARAMETERS

write(*,*),"    Reading parameters "
open(unit=10,file='input.dat',status='old')
 read(10,*),m,nb,am,ncyc,PF,PT,T
 NP=ncyc/PF
 write(2,*), NP
 Beta=1/(T*8.617343d-5 )
 read(10,*),ntypes,nid,njumps
 allocate(barr(njumps,nid,nid,ntypes,ntypes,ntypes,ntypes))
 allocate(freqs(njumps,nid,nid,ntypes,ntypes,ntypes,ntypes))
 allocate(twobody(njumps,nid,nid,ntypes,ntypes))
 allocate(amoves(njumps,nid,nid,ntypes,ntypes))
 allocate(onebody(nid,ntypes))
 onebody=0.0
 twobody=0.0
 amoves=-1.0
 barr=0.0
 freqs=0.0
 do i=1,m
  read(10,*),jl,id1,id2,t1,t2,E
  twobody(jl,id1,id2,t1,t2)=E
  twobody(jl,id2,id1,t2,t1)=E
 enddo

 do i=1,nb
  read(10,*) id1,t1,E
  onebody(id1,t1)=E
 enddo

 do i=1,am
  read(10,*),jl,id1,id2,t1,t2,nt1,nt2,v0,E
  amoves(jl,id1,id2,t1,t2)=1.0
  amoves(jl,id2,id1,t2,t1)=1.0
  if (E > 0) then
   barr(jl,id1,id2,t1,t2,nt1,nt2)=E
   freqs(jl,id1,id2,t1,t2,nt1,nt2)=v0
  endif

 enddo

write(*,*), "    Complete"
write(*,*), "    ---------------------------------------------"



 do i=1,n                 !Update all jumps
  c1=0
  sites(i)%njump=c1
  do j=1,sites(i)%nons
    do nt1=1,ntypes
     do nt2=1,ntypes

    in=sites(i)%ns(j)
    jl=sites(i)%jt(j)
    id1=sites(i)%id
    id2=sites(in)%id
    t1=sites(i)%species
    t2=sites(in)%species

   if ( barr(jl,id1,id2,t1,t2,nt1,nt2) > 0.0 ) then
   c1=c1+1
   sites(i)%jw(c1)=sites(i)%ns(j)
   sites(i)%njump=c1
   sites(i)%Ea0(c1)=barr(jl,id1,id2,t1,t2,nt1,nt2)
   sites(i)%p(c1)=0.0
   sites(i)%freq(c1)=freqs(jl,id1,id2,t1,t2,nt1,nt2)
   sites(i)%eff(c1,1)=t1
   sites(i)%eff(c1,2)=t2
   sites(i)%eff(c1,3)=nt1
   sites(i)%eff(c1,4)=nt2

   endif

     enddo
    enddo
  enddo
 enddo                      !/update all jumps



! #E# READ IN BARR

! #S# CALCULATE TOTAL ENERGY
Etot=0.0
do i=1,n
 Etot=Etot+2*onebody(sites(i)%id,sites(i)%species) ! to balance the doulbe counting in the 2-body
 do j=1,sites(i)%nons
Etot=Etot+twobody( sites(i)%jt(j), sites(i)%id,  &
& sites( sites(i)%ns(j))%id, sites(i)%species, &
& sites( sites(i)%ns(j) )%species )
 enddo
enddo
Etot=Etot*0.5
write (*,*) "Initial total energy: ",Etot
! #E# /CALCULATE TOTAL ENERGY
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!








!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! #S# MAIN
counter1=0
do i=1,ncyc


!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! UPDATING ALL PROBABILITIES
! LOOP OF ALL POSSIBLE HOPS
do ii=1,n
 hoppf=ii
 do jj=1,sites(hoppf)%njump
  hoppt=sites(hoppf)%jw(jj)

! DETERMINE THE FINAL STATES
fsites=sites
fsites(hoppf)%species=sites(hoppf)%eff(jj,3)
fsites(hoppt)%species=sites(hoppf)%eff(jj,4)
! /DETERMINE THE FINAL STATES

! DETERMINE ENERGY OF A AND B IN INITIAL STATE
Einit=0
do j=1,sites(hoppf)%nons
Einit=Einit+twobody( sites(hoppf)%jt(j), sites(hoppf)%id,  &
& sites( sites(hoppf)%ns(j))%id, sites(hoppf)%species, &
& sites( sites(hoppf)%ns(j) )%species )
enddo

do j=1,sites(hoppt)%nons
Einit=Einit+twobody( sites(hoppt)%jt(j), sites(hoppt)%id, &
& sites( sites(hoppt)%ns(j))%id, sites(hoppt)%species, &
& sites(sites(hoppt)%ns(j))%species )
enddo

Einit=Einit+onebody( sites(hoppf)%id, sites(hoppf)%species)
Einit=Einit+onebody( sites(hoppt)%id, sites(hoppt)%species)
! /DETERMINE ENERGY OF A AND B IN INITIAL STATE

! DETERMINE ENERGY OF A AND B IN FINAL STATE
Efinal=0
do j=1,fsites(hoppf)%nons
Efinal=Efinal+twobody( fsites(hoppf)%jt(j), fsites(hoppf)%id,  &
& fsites( fsites(hoppf)%ns(j))%id, fsites(hoppf)%species, &
& fsites( fsites(hoppf)%ns(j))%species )
enddo

do j=1,fsites(hoppt)%nons
Efinal=Efinal+twobody( fsites(hoppt)%jt(j), fsites(hoppt)%id,  &
& fsites( fsites(hoppt)%ns(j))%id, fsites(hoppt)%species, &
& fsites( fsites(hoppt)%ns(j))%species )
enddo

Efinal=Efinal+onebody( fsites(hoppt)%id, fsites(hoppt)%species)
Efinal=Efinal+onebody( fsites(hoppf)%id, fsites(hoppf)%species)


! /DETERMINE ENERGY OF A AND B IN FINAL STATE
dE=Efinal-Einit
Ea=de2ea(dE, sites(hoppf)%Ea0(jj) )
sites(ii)%p(jj) = sites(ii)%freq(jj)*exp(- Ea*Beta )
enddo
enddo
!/LOOP OVER ALL POSSIBLE HOPS
!/UPDATING PROBABILITIES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Calculate tot prob and decide jump
totp=0.0
do ii=1,n
  do j=1,sites(ii)%njump
    totp=totp+sites(ii)%p(j)
  enddo
enddo

r=ran()*totp
dtime=-log(ran())/totp
totp=0.0


do ii=1,n
  do j=1,sites(ii)%njump
    totp=totp+sites(ii)%p(j)
    if (totp > r) then
      go to 10
    endif
  enddo
enddo

10 time=time+dtime
!/Calculate tot prob and decide jump
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!UPDATE GEOMETRY
sites(ii)%species=sites(ii)%eff(j,3)
sites( sites(ii)%jw(j) )%species=sites(ii)%eff(j,4)

i1=ii
i2=sites(ii)%jw(j)

  do l=1,sites(i1)%nons ! First involed and its ns
  ii=sites(i1)%ns(l)
  c1=0
  sites(ii)%njump=c1
  do j=1,sites(ii)%nons
    do nt1=1,ntypes
     do nt2=1,ntypes

    in=sites(ii)%ns(j)
    jl=sites(ii)%jt(j)
    id1=sites(ii)%id
    id2=sites(in)%id
    t1=sites(ii)%species
    t2=sites(in)%species

   if ( barr(jl,id1,id2,t1,t2,nt1,nt2) > 0.0 ) then
   c1=c1+1
   sites(ii)%jw(c1)=sites(ii)%ns(j)
   sites(ii)%njump=c1
   sites(ii)%Ea0(c1)=barr(jl,id1,id2,t1,t2,nt1,nt2)
   sites(ii)%freq(c1)=freqs(jl,id1,id2,t1,t2,nt1,nt2)
   sites(ii)%p(c1)=0.0
   sites(ii)%eff(c1,1)=t1
   sites(ii)%eff(c1,2)=t2
   sites(ii)%eff(c1,3)=nt1
   sites(ii)%eff(c1,4)=nt2

   endif

     enddo
    enddo
  enddo
  enddo


  do l=1,sites(i2)%nons ! Second involed and ns
  ii=sites(i2)%ns(l)
  c1=0
  sites(ii)%njump=c1
  do j=1,sites(ii)%nons
    do nt1=1,ntypes
     do nt2=1,ntypes

    in=sites(ii)%ns(j)
    jl=sites(ii)%jt(j)
    id1=sites(ii)%id
    id2=sites(in)%id
    t1=sites(ii)%species
    t2=sites(in)%species

   if ( barr(jl,id1,id2,t1,t2,nt1,nt2) > 0.0 ) then
   c1=c1+1
   sites(ii)%jw(c1)=sites(ii)%ns(j)
   sites(ii)%njump=c1
   sites(ii)%Ea0(c1)=barr(jl,id1,id2,t1,t2,nt1,nt2)
   sites(ii)%freq(c1)=freqs(jl,id1,id2,t1,t2,nt1,nt2)
   sites(ii)%p(c1)=0.0
   sites(ii)%eff(c1,1)=t1
   sites(ii)%eff(c1,2)=t2
   sites(ii)%eff(c1,3)=nt1
   sites(ii)%eff(c1,4)=nt2

   endif

     enddo
    enddo
  enddo
  enddo
!/UPDATE GEOMETRY
!!!!!!!!!!!!!!!!!!


! #S# PRINT
if ( MOD(i,PF) .EQ. 0 ) then
 write(2,*),n
 write(*,*),time,dtime
 do pi=1,n
 write(2,*),sites(pi)%species,sites(pi)%x,sites(pi)%y,sites(pi)%z
 enddo
endif
! #E# PRINT



enddo
!/END OF MAIN
!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!


! #S# PRINT geom.out
write(3,*), n
 do i=1,n
write(3,*),sites(i)%x, &
sites(i)%y, &
sites(i)%z, &
sites(i)%id, &
sites(i)%species, &
sites(i)%nons, &
(sites(i)%ns(j), j=1,sites(i)%nons), &
(sites(i)%jt(j), j=1,sites(i)%nons)
enddo



! #E# PRINT geom.out



close(2)
close(3)
end program MC

! #E#


function de2ea(E,Ea0)
  double precision, intent(in) :: E,Ea0
  double precision :: de2ea

   if(E <= 0) then
   de2ea=Ea0+E
    if(de2ea < 0.0) then
    de2ea=0.0
    endif
   endif

   if(E > 0) then
   de2ea=E
    if(de2ea < Ea0) then
    de2ea=Ea0
    endif
   endif


end function





!     - -         - -   - -   - -         - -         - -   - -   - -
!    |     |   |   | | |   | |   | |   |   |     |   |   | |     |
!     - -           -   - -                                 - -   - -
!        | |   |   | | |  \  |   | |   |   |     |   |   | |         |
!     - -   - -   - -         - -   - -                     - -   - -


subroutine print_out(sites,n)
type site
  double precision :: x,y,z
  integer          :: id               ! type of site
  integer          :: species          ! type of species (vacancy , normal oxygen, ad-atom etc.)
  integer          :: nons             ! no of neigbors
  integer          :: njump            ! no of possible jumps
  integer, dimension(99) :: ns ! nieghbohring sites
  integer, dimension(99) :: jt ! jumptype to ns
  double precision, dimension(99)   :: p ! probability of the jump
  integer, dimension(99)            :: jw ! jump where
  integer, dimension(99,4)          :: eff  ! effect of the jump (njump),(site1,site2,newspec1,newspec2)
endtype
type(site), dimension(n), intent(in) :: sites
integer :: i,n


do i=1,n
write(*,*),sites(i)%species,sites(i)%x,sites(i)%y,sites(i)%z
enddo

end subroutine print_out
