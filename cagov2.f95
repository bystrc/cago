program cago
  !----------------------------------------
  use pdbx
  !----------------------------------------
  implicit none
  !---------------------------------------
  include "mpif.h"
  integer,parameter::traj_unit=98 !unique file handle for trajectory file
  integer,parameter::sample_unit=16 !unique file handle for energy log file
  integer,parameter::ddihed_unit=457 ! unique file handle for dihedral monitor log file
  integer,parameter::nlocal=4
  integer,parameter::unfold_step=30000 !predetermined no. of iterations
                                       !for unfolding
  !----------------------------------------
  real,parameter::pi=3.14159265
  real,parameter::vdwr=4.0 !vanderWaal's radius
  real,parameter::vdwr12=vdwr**12 !1st term in Lennard-Jones 6-12
  real,parameter::deg2rad=pi/180.
  real,parameter::rad2deg=180/pi
  real,parameter::dt=0.05
  real,parameter::tcoup=0.3
  real,parameter::unfold_temp=3.0 !preset unfolding temperature (canonical)
  real,parameter::friction=0.1 !to slow down molecule after hit with
                               !random force vector
  !----------------------------------------
  integer,parameter::k_max_consec_rmsd=200 !if the rmsd is consistent for this many frames then the program will stop
  integer,parameter::k_max_consec_q=200 !if the q is greater than 0.99 for this many frames then the program will stop
  real,parameter::k_acptble_rmsd_range=1.3 !this is the acceptable tolerance for difference in rmsd from the comparison rmsd
  real,parameter::k_bond=100 !converesion factor for bond stretch/contraction to energy     
  real,parameter::k_anglebend=20 !conversion factor for angle bending (theta) to
                                 !energy
  real ::k_dih1  !conversion factors dih1 and 3 for dihedral to energy
  real ::k_dih3
  real :: expo !slope steepness for dihedral barrier
  real :: min_q  ! minimum q value for monitoring unfoldedness
  real,parameter::k_contact=1.0
  real,parameter::k_ljrep=0.1
  !real,parameter::contact_cut=14.0
  !real,parameter::calpha_tol=8.0
  real,parameter::allatm_tol=6.5
  !----------------------------------------
  type(residue),dimension(:),pointer::protein
  !----------------------------------------
  
integer::codelen
character(len=1)::chain_id
  character(len=4)::code
  character(len=100)::nrg
  character(len=6)::pro_num
  character(len=100)::mpdb
  character(len=100)::pdb
  character(len=100)::aline
  character(len=4)::test
  character(len=100)::test2
  integer::n
  integer::jarg
  integer::ierr
  integer::nres
  integer::ndsbonds !number of disulfide bonds
  integer::go_s !starting residue for GO
  integer::go_e !ending residue for GO
  integer::istep
  integer::nstep
  integer::sstep
  integer::reclen
  integer::num_consec_rmsd, num_consec_q
  integer::ios
  integer,dimension(3000)::iseq
  integer,dimension(:),allocatable,target::seed
  real::sys_temp
  real::fold_temp
  real::rep_rmsd !rmsd that gets checked against future rmsds for consistency (never printed)
  real::q
  real::rmsd
  real::curr_temp
  real::t1
  real::t2
  real::bond_e
  real::total_e
  real::anglebend_e
  real::dihedral_e
  real::contact_e
  real::ljrep_e
  real::dsbond_e !disulfide bond energy
  real,dimension(2)::tarry
  real,dimension(:),pointer::target_anglebend
  real,dimension(:),pointer::target_dihedral
  real,dimension(:),pointer::angd !stores the current SS dihedral
  real,dimension(:),pointer::angp !stores the previous SS dihedral calculations
  real,dimension(:),pointer::rotd !stores the direction of rotation from angp &
! to angd
  real,dimension(:),pointer::sumr !stores overall rotation at the end of the run
  !real,dimension(:),allocatable::prot
  real,dimension(:,:),pointer::native
  real,dimension(:,:),pointer::model
  integer,dimension(:,:),pointer::dsbonds
  real,dimension(:,:),pointer::net_force
  real,dimension(:,:),pointer::anglebend_f
  real,dimension(:,:),pointer::dihedral_f
  real,dimension(:,:),pointer::contact_f
  real,dimension(:,:),pointer::cmap
  real,dimension(:,:),pointer::bond_f
  real,dimension(:,:),pointer::velocity
  real,dimension(:,:),pointer::random_f
  real,dimension(:,:),pointer::ljrep_f
  real,dimension(:,:),pointer::acceleration
  real,dimension(:,:),pointer::dsbond_f
  !------------------------------------------------------------------------------------
  integer::countl
  integer::ires
  integer::seed_1
  integer::seed_2
  integer::seed_3
  integer::seed_4
  !integer::seed_5
  !integer::seed_6
  !integer::seed_7
  !integer::seed_8
  logical::unfolded
  logical::currently_knotted
  integer::c_p
  integer::t_p
  integer::i_p
  integer::length
  integer::frame
  integer :: run_id
  character(len=8)::num
  integer,dimension(MPI_STATUS_SIZE)::stat
  logical :: movieout=.true.
  !------------------------------------------------------------------------------------
    jarg = command_argument_count()
  if(jarg < 11) then
    write(*,*)
    write(*,*)'usage: xcago {native pdb} {chain id} {protein code} {temp}'// &
     " {nstep} {sample step} {dihedral_k} {dihedral_expo} {run ID} {Start residue} {End residue} "
    write(*,*)
    stop
  endif
  call flush(6)
  !------------------------------------------------------------------------------------
  call getarg(1,pdb)
  call getarg(2,chain_id)
  call getarg(3,code)
  call getarg(4,aline)
  !------------------------------------------------------------------------------------
  read(aline,*,iostat=ierr)fold_temp
  if(ierr /= 0) then
    write(*,*) "error xcago :: fold temp is not a number! ",fold_temp
  endif
  call getarg(5,aline)
  read(aline,*,iostat=ierr)nstep
  if(ierr /= 0) then
    write(*,*) "error xcago :: nsteps is not a number! ",nstep
  endif
  call getarg(6,aline)
  read(aline,*,iostat=ierr)sstep
  if(ierr /= 0) then
    write(*,*) "error xcago :: sample_step is not a number! ",sstep
  endif
  call getarg(7,aline)
  read(aline,*,iostat=ierr)k_dih1
  if(ierr /= 0) then
    write(*,*) "dihedral_k is not a number"
  endif
  call getarg(8,aline)
  read(aline,*,iostat=ierr) expo
  if(ierr /= 0) then
    write(*,*) "dihedral exponent is not a number"
  endif
  call getarg(9,aline)
  read(aline,*,iostat=ierr)run_id
  if(ierr /= 0) then
    write(*,*) "Run ID is not a number"
  endif 
  call getarg(10,aline)
  read(aline,*,iostat=ierr)go_s
  if(ierr /= 0) then
    write(*,*) "error xcago :: residue number is not an integer! ",sstep
  endif
  call getarg(11,aline)
  read(aline,*,iostat=ierr)go_e
  if(ierr /= 0) then
    write(*,*) "error xcago :: residue number is not an integer! ",sstep
  endif
  call flush(6)
  !if(jarg > 7) then
  !  call getarg(8,aline)
  !  read(aline,*,iostat=ierr)friction
  !  if(ierr /= 0) then
  !    write(*,*) "error xcago :: friction is not a number! ",friction
  !  endif
  !else
  !  friction=0.1
  !endif
  call pdbx_read_protein(pdb,chain_id,protein,dsbonds,ndsbonds,iseq,native,nres)
  call get_cmap(protein,nres,cmap)
  !------------------------------------------------------------------------------------
  call MPI_INIT(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD,c_p,ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,t_p,ierr)
  
  run_id = run_id*10000 + c_p+fold_temp*1000
  write(pro_num,'(I6.6)')run_id

  nrg=TRIM(code)//pro_num// ".nrg"
  mpdb=TRIM(code)//pro_num//".mpdb"
  open(sample_unit,file=nrg,iostat=ierr)
  if(ierr /= 0) stop 'error xcago :: cannot open sample file'
  open(ddihed_unit,file=TRIM(code)//pro_num//".log", iostat=ierr) !attempt to open dihedral log file
  if(ierr/=0) stop 'error xcago :: cannot open dihedral log file'  
  !------------------------------------------------------------------------------------
  if (movieout) then
    if(traj_unit /= 99) then
      open(traj_unit,file=mpdb,status='replace',iostat=ierr)
      if(ierr /= 0) stop 'Cannot open trajectory file!'
    else
      reclen=12*nres+8
      open(traj_unit,file=mpdb,access='direct',form='unformatted',recl=reclen,iostat=ierr,status='replace')
      if(ierr /= 0) stop 'Cannot open trajectory file!'
    endif
  endif
 
  !------------------------------------------------------------------------------------
  call random_seed(size=n)
  allocate(seed(n))
  if(c_p == 0) then
  do i_p=1,t_p-1
   call init_ranseed(i_p)
   call random_seed(get=seed) 
    call MPI_SEND(seed(1),n,MPI_INTEGER,i_p,1,MPI_COMM_WORLD,ierr)
  enddo
  endif
  if(c_p == 0) then
   call init_ranseed(0)
   call random_seed(get=seed)
  endif
  if(c_p /= 0) then
    call MPI_RECV(seed(1),n,MPI_INTEGER,0,1,MPI_COMM_WORLD,stat,ierr)
    call random_seed(put=seed)
  endif
  call flush(6)
  !------------------------------------------------------------------------------------
  write(sample_unit,*)"---------------------------------------------------------------------"
  write(sample_unit,'(1x,a)')"Performing simulation with Reduced Units (Lennard-Jones)... mass=epilson=sigma=k_boltz=1"
  write(sample_unit,'(1x,a,f10.3)'  )"Fold Temperature  ",fold_temp
  write(sample_unit,'(1x,a,f10.3)'  )"Friction            ",friction
  write(sample_unit,'(1x,a,f10.3)'  )"Timestep            ",dt
  write(sample_unit,'(1x,a,i10)'    )"Nsteps              ",nstep
  write(sample_unit,'(1x,a,i10)'    )"Sample Step         ",sstep
  write(sample_unit,'(1x,a,a,1x,a)' )"Native{chain}       ",trim(pdb),chain_id
  write(sample_unit,'(1x,a,i10)'    )"Nres                ",nres
  write(sample_unit,'(1x,a,8i12)'   )"Seed                ",seed
  write(sample_unit,*)"---------------------------------------------------------------------"
  !------------------------------------------------------------------------------------
  !if((go_s/=0 .or. go_e/=0) .and. go_e>go_s) then
  !nres=go_e-go_s+1
  !endif
  if(.not. associated(ljrep_f)) allocate(       ljrep_f(3,nres))
  if(.not. associated(model))allocate(         model(3,nres))
  if(.not. associated(net_force))allocate(     net_force(3,nres))
  if(.not. associated(bond_f))allocate(        bond_f(3,nres))
  if(.not. associated(random_f))allocate(      random_f(3,nres))
  if(.not. associated(velocity))allocate(      velocity(3,nres))
  if(.not. associated(anglebend_f))allocate(   anglebend_f(3,nres))
  if(.not. associated(dihedral_f))allocate(    dihedral_f(3,nres))
  if(.not. associated(contact_f))allocate(     contact_f(3,nres))
  if(.not. associated(acceleration))allocate(  acceleration(3,nres))
  if(.not. associated(target_dihedral))allocate( target_dihedral(nres))
  if(.not. associated(target_anglebend))allocate(target_anglebend(nres))
  if(.not. associated(dsbond_f))allocate(      dsbond_f(3,nres))
  write(*,*) "all arrays allocated"
  write(*,*) "Disulfide bonds:", dsbonds
  call flush(6)
  !------------------------------------------------------------------------------------
  model=native
  !do ires=1,3
  ! do jres=0,nres-1
  !model(ires,jres+1)=prot(i_res,j_res+go_s)
  !native=model
  sys_temp=unfold_temp
  !if(ndsbonds>0) then
    allocate(angd(2))
    write (*,*)'Initializing rotd'
    allocate(rotd(2))
    rotd=0. !initialize rotation metric
    write (*,*)'Initialized rotd. initializing sumr '
    allocate(sumr(2))
    sumr=0  !initialize overall rotation metric
  !endif
  !------------------------------------------------------------------------------------
  !call rand_chain(model,nres) !comment out to start at native state 
  call update_anglebend(target_anglebend,native,nres)
  !write (*,*) 'updated angle bend'
  call update_dihedral(target_dihedral,native,nres)
  !write(*,*) 'dihedral'
  call set_target_dihedral(target_dihedral,nres)
  !write(*,*)'target dihedral set'
  if (movieout) call write_traj(traj_unit,model,nres,istep,iseq,sys_temp)
  !write(*,*)'trajectory 1 written'
  call bond_nrg(native,nres,bond_e)
  !write(*,*)'initial bond energy calculated'
  call dsbond_nrg(dsbonds,native,model,nres,ndsbonds,dsbond_e)
  !write(*,*)'Initial disulfide energy calculated ndsbonds=',ndsbonds
  if(associated(angd)) deallocate(angd)
  allocate(angd(2))
  !if(ndsbonds>0) then
    call sdhmon(angd,native,dsbonds)
    if (associated(angp)) deallocate(angp)
    allocate(angp(2))
    angp=angd !initialize the previous SS dihedral array--one-time
    call sdhdir(angd,angp,rotd,sumr)
  !endif
  call anglebend_nrg(native,anglebend_e,nres,target_anglebend)
  !write(*,*)'Angle bend calculated'
  call dihedral_nrg(native,dihedral_e,nres,target_dihedral)
  !write(*,*)'dihedral energy'
  !call lj_repulsion_nrg(native,ljrep_e,nres)
  call contact_nrg(native,contact_e,cmap,nres)
  !write(*,*)'contact nrg'
  call get_rmsd(native,native,rmsd,nres)
  !write(*,*)'rmsd'
  call get_Q(native,cmap,nres,q)
  !write(*,*)'Q value obtained, starting MD'
  !---------------------------------------------------------------
  total_e = bond_e + anglebend_e + dihedral_e + contact_e + dsbond_e !+ ljrep_e
  !---------------------------------------------------------------
  !write(sample_unit,'(i10,8f10.2,2f6.2)')0,total_e,bond_e,dsbond_e,anglebend_e,dihedral_e,contact_e,sys_temp,q,rmsd
  write(sample_unit,*)0,total_e,bond_e,dsbond_e,anglebend_e,dihedral_e,contact_e,sys_temp,q,rmsd
  !write(*,'(i10,6f10.2,2f6.2)')0,total_e,bond_e,anglebend_e,dihedral_e,contact_e,sys_temp,q,rmsd
  !------------------------------------------------------------------------------------
  call bond_force(model,nres,bond_f)
  !write (*,*) '1'
  call dsbond_force(model,native,dsbonds,ndsbonds,nres,dsbond_f)
  !write (*,*) '2'
  call anglebend_force(model,anglebend_f,nres,target_anglebend)
  !write (*,*) '3'  
  call dihedral_force(model,dihedral_f,nres,target_dihedral)
  !write (*,*) '4'
  !call lj_repulsion_force(model,ljrep_f,nres)
  call contact_force(model,cmap,contact_f,nres)
  !write (*,*) '5'
  call random_force(model,random_f,sys_temp,friction,dt,nres)
  !write (*,*) '6'
  call initialize_velocity(model,velocity,nres,sys_temp)
  !------------------------------------------------------------------------------------
  net_force = bond_f + contact_f + dihedral_f + anglebend_f + dsbond_f !+ ljrep_f
  !------------------------------------------------------------------------------------
  acceleration = net_force - friction*velocity + random_f 
  !---------------------------------------------------------------
  call etime(tarry,t1)
  !------------------------------------------------------------------------------------
  rep_rmsd = -1 !set to impossible value for initial initialization
  num_consec_rmsd = 0
  num_consec_q = 0
  unfolded= .false.
  min_q = 1.00
  do istep=1,nstep
      ! if(num_consec_rmsd >= k_max_consec_rmsd .and. rep_rmsd < 3) exit !checks to see if proteins structure is more or less stable by looking for changing rmsd
      if(num_consec_q >= k_max_consec_q ) exit !checks to see if protein is folded
      !---------------------------------------------------------------
  
!if (istep > unfold_step .and. mod(istep,100)==0) sys_temp=sys_temp*0.9


    if(istep > unfold_step) sys_temp=0.133
    !if(istep > unfold_step) sys_temp=fold_temp
      !---------------------------------------------------------------
      call velocity_verlet_1(acceleration,model,velocity,nres,dt) !change positions and velocities
      !---------------------------------------------------------------
      call bond_force(model,nres,bond_f)
      call dsbond_force(model,native,dsbonds,ndsbonds,nres,dsbond_f)
      call anglebend_force(model,anglebend_f,nres,target_anglebend)
      call flush(6)
      call dihedral_force(model,dihedral_f,nres,target_dihedral)
      !call lj_repulsion_force(model,ljrep_f,nres)
      call contact_force(model,cmap,contact_f,nres)
      call random_force(model,random_f,sys_temp,friction,dt,nres)
      !---------------------------------------------------------------
      net_force = bond_f + contact_f + anglebend_f + dihedral_f + dsbond_f !+ ljrep_f
      !---------------------------------------------------------------
      acceleration = net_force - friction*velocity + random_f 
      !---------------------------------------------------------------
      call velocity_verlet_2(acceleration,velocity,nres,dt) !modify velocities
      !---------------------------------------------------------------
      call berendsen_tcoup(sys_temp,curr_temp,velocity,nres)
      !---------------------------------------------------------------
      call sdhmon(angd,model,dsbonds)
           !write (*,*) angp(:),'----',angd(:)
           call sdhdir(angp,angd,rotd,sumr)
           if(.not. associated(rotd)) then
              write(*,*)'Error! rotd not allocated'
           endif
           if(.not. associated(sumr)) then
              write(*,*)'Error! rotd not allocated'
           endif
      write (ddihed_unit,*) rotd,'---',sumr !write out dihedral data to a file 
      call flush(ddihed_unit)
      if(mod(istep,sstep)==0) then
        call bond_nrg(model,nres,bond_e)
        call dsbond_nrg(dsbonds,native,model,nres,ndsbonds,dsbond_e)
        !if(ndsbonds > 0) then
          ! call sdhmon(angd,model,dsbonds)
          ! !write (*,*) angp(:),'----',angd(:)
          ! call sdhdir(angp,angd,rotd,sumr)
          ! if(.not. associated(rotd)) then
          !    write(*,*)'Error! rotd not allocated'
          ! endif
          ! if(.not. associated(sumr)) then
          !    write(*,*)'Error! rotd not allocated'
          ! endif
        !endif
        call anglebend_nrg(model,anglebend_e,nres,target_anglebend)
        call dihedral_nrg(model,dihedral_e,nres,target_dihedral)
!        !call lj_repulsion_nrg(model,ljrep_e,nres)
        call contact_nrg(model,contact_e,cmap,nres)
        call get_rmsd(native,model,rmsd,nres)
        call get_Q(model,cmap,nres,q)
        call remove_net_momentum(model,velocity,nres)
        if (movieout) call write_traj(traj_unit,model,nres,istep,iseq,sys_temp)
        !---------------------------------------------------------------
        total_e = bond_e + anglebend_e + dihedral_e + contact_e + dsbond_e !+ ljrep_e
        if(total_e >= 10000)then
          write(*,*) "process ",c_p,"(",trim(mpdb),") blew up: energy:",total_e
          call MPI_FINALIZE(ierr)
          stop
        endif
        !---------------------------------------------------------------
        !write(sample_unit,'(i10,7f10.2,2f6.2)')istep,total_e,bond_e,anglebend_e,dihedral_e,contact_e,ljrep_e,sys_temp,q,rmsd
        !write(sample_unit,'(i10,7f10.2,2f6.2,2f6.6,a,2f6.6,a,2f3.5)')&
        !     istep,total_e,bond_e,dsbond_e,anglebend_e,dihedral_e &
        !     ,contact_e,sys_temp,q,rmsd,angd(:),'-',rotd(:),'-',sumr(:)
        !write(sample_unit,*) istep,total_e,bond_e,dsbond_e,anglebend_e,dihedral_e &
        !     ,contact_e,sys_temp,q,rmsd,angd(:),'-',rotd(:),'-',sumr(:) 
        write(sample_unit,*) istep,total_e,bond_e,dsbond_e,anglebend_e,dihedral_e &
             ,contact_e,sys_temp,q,rmsd 
        !write(*,'(i10,6f10.2,2f6.2)')istep,total_e,bond_e,anglebend_e,dihedral_e,contact_e,sys_temp,q,rmsd
        !---------------------------------------------------------------
        if (q > 0.99) then
          num_consec_q = num_consec_q + 1
        else
          num_consec_q = 0
          if (q < min_q) min_q = q
        endif
      endif
      call flush(sample_unit)
  enddo
  write(*,*) "FINISHED on process ",c_p," Minimum q=",min_q," Temperature=",sys_temp," Unfolding time=",istep-unfold_step
  call etime(tarry,t2)
  write(sample_unit,*)(t2-t1)/60.0
  !------------------------------------------------------------------------------------
  if(associated(model))               deallocate(model)
  if(associated(net_force))           deallocate(net_force)
  if(associated(target_dihedral))     deallocate(target_dihedral)
  if(associated(target_anglebend))    deallocate(target_anglebend)
  if(associated(bond_f))              deallocate(bond_f)
  if(associated(random_f))            deallocate(random_f)
  if(associated(anglebend_f))         deallocate(anglebend_f)
  if(associated(dihedral_f))          deallocate(dihedral_f)
  if(associated(acceleration))        deallocate(acceleration)
  if(associated(ljrep_f))             deallocate(ljrep_f)
  if(associated(dsbonds))             deallocate(dsbonds)
  if(associated(dsbond_f))            deallocate(dsbond_f)
  do countl=1,ndsbonds
     !if(sumr(countl)>0) then
       write(sample_unit,*) 'The overall rotation for disulfide bond '&
       ,countl,' was ',sumr(countl),'degrees'
     !endif
     !if(sumr(countl)<0) then
     !  write(sample_unit,'(1x,a,i2,a,f6.6,a)') 'The overall rotation for disulfide bond ' &
     !  ,countl,' was ',abs(sumr(countl)),'degrees to the left'
     !endif
     !if(sumr(countl)==0) then
     !  write(sample_unit,'(1x,a,i2)'),'No rotation observed for disulfide bond ',countl
     !endif
  enddo
  call flush(sample_unit)
  if (movieout) close(traj_unit)
  close(sample_unit)
  close(ddihed_unit)
  !------------------------------------------------------------------------------------
 call MPI_FINALIZE(ierr)
 contains
  logical function isnan(x)
    implicit none
    real,intent(in) :: x
    isnan = .false.
    if( x /= x) isnan = .true.
  end function isnan

!----------------------------------------------------------------------------------------------------------------------------
!  subroutine  lj_repulsion_nrg(xyz,nrg,nres)
!   !---------------------------------------
!   implicit none
!   !---------------------------------------
!   integer,intent(in)::nres
!   real,intent(out)::nrg
!   real,dimension(:,:),intent(in)::xyz
!   !---------------------------------------
!   integer::ires
!   integer::jres
!   real::x
!   real::dij
!   real::dijsq
!   real::dij3
!   real::dij6
!   real::dij12
!   real,dimension(3)::rij
!   !---------------------------------------
!   nrg=0
!   do ires = 1, nres-nlocal
!     do jres = ires + nlocal, nres
!         rij = xyz(:,ires) - xyz(:,jres)
!         dijsq = sum(rij*rij)
!         dij = sqrt(dijsq); dij3=dijsq*dij; dij6=dij3*dij3; dij12=dij6*dij6 
!         x = k_ljrep*((vdwr12/dij12)) 
!         nrg = nrg + x
!     enddo
!   enddo
!   !---------------------------------------
!  end subroutine  lj_repulsion_nrg
!!----------------------------------------------------------------------------------------------------------------------------
!  subroutine  lj_repulsion_force(xyz,lj_f,nres)
!   !---------------------------------------
!   implicit none
!   !---------------------------------------
!   integer,intent(in)::nres
!   real,dimension(:,:),intent(in)::xyz
!   real,dimension(:,:),intent(out)::lj_f
!   !---------------------------------------
!   integer::ires
!   integer::jres
!   real::x
!   real::dij
!   real::dijsq
!   real::dij3
!   real::dij6
!   real::dij12
!   real,dimension(3)::rij
!   !---------------------------------------
!   lj_f=0
!   do ires = 1, nres-nlocal
!     do jres = ires + nlocal, nres
!         rij = xyz(:,ires) - xyz(:,jres)
!         dijsq = sum(rij*rij)
!         dij = sqrt(dijsq); dij3=dijsq*dij; dij6=dij3*dij3; dij12=dij6*dij6
!         x = k_ljrep*12*(vdwr12/(dij12*dij)) 
!         lj_f(:,ires) = lj_f(:,ires) + x*rij/dij
!         lj_f(:,jres) = lj_f(:,jres) - x*rij/dij
!     enddo
!   enddo
!   !---------------------------------------
!  end subroutine  lj_repulsion_force
!----------------------------------------------------------------------------------------------------------------------------
  subroutine contact_nrg(xyz,nrg,cmap,nres)
   implicit none
   integer,intent(in)::nres
   real,intent(out)::nrg
   real,dimension(:,:),intent(inout)::xyz
   real,dimension(:,:),pointer::cmap
   real::x,dij,dijsq,dij3,dij6,dij10,dij12,sig12
   real,dimension(3)::rij
   integer::ires,jres

   nrg=0.
   do ires = 1, nres-nlocal
     do jres = ires + nlocal, nres
       if(cmap(ires,jres)/=0)then
         rij = xyz(:,ires) - xyz(:,jres)
         dijsq = sum(rij*rij)
         dij = sqrt(dijsq); dij3=dijsq*dij; dij6=dij3*dij3; dij12=dij6*dij6; dij10=dijsq*dijsq*dij6
         sig12=cmap(jres,ires)*cmap(ires,jres)*cmap(ires,jres)
         x = k_contact*(5*(sig12/dij12)-6*(cmap(jres,ires)/dij10))
         nrg = nrg + x
       else  
         rij = xyz(:,ires) - xyz(:,jres)
         dijsq = sum(rij*rij)
         dij= sqrt(dijsq); dij3=dijsq*dij; dij6=dij3*dij3; dij12=dij6*dij6
         x = k_ljrep*vdwr12/dij12
         nrg = nrg + x
       endif
     enddo
   enddo
  end subroutine
!----------------------------------------------------------------------------------------------------------------------------
  subroutine contact_force(xyz,cmap,contactf,nres)
   implicit none
   integer,intent(in)::nres
   real,dimension(:,:),intent(inout)::xyz
   real,dimension(:,:),pointer::cmap,contactf
   real::x,dij,dijsq,dij3,dij6,dij10,dij12,sig12
   real,dimension(3)::rij
   integer::ires,jres

   contactf=0
   do ires = 1, nres-nlocal
     do jres = ires + nlocal, nres
       !if(ires/=go_s .and. ires/=go_e) then
       if(cmap(ires,jres)/=0)then
         rij = xyz(:,ires) - xyz(:,jres)
         dijsq = sum(rij*rij)
         dij = sqrt(dijsq); dij3=dijsq*dij; dij6=dij3*dij3; dij12=dij6*dij6; dij10=dijsq*dijsq*dij6
         sig12=cmap(jres,ires)*cmap(ires,jres)*cmap(ires,jres)
         x = k_contact*60*(sig12/(dij*dij12)-cmap(jres,ires)/(dij10*dij))
         contactf(:,ires) = contactf(:,ires) + x*rij/dij    
         contactf(:,jres) = contactf(:,jres) - x*rij/dij 
       else  
         rij = xyz(:,ires) - xyz(:,jres)
         dijsq = sum(rij*rij)
         dij= sqrt(dijsq); dij3=dijsq*dij; dij6=dij3*dij3; dij12=dij6*dij6
         x = k_ljrep*12*vdwr12/(dij12*dij)
         contactf(:,ires) = contactf(:,ires) + x*rij/dij
         contactf(:,jres) = contactf(:,jres) - x*rij/dij
       endif
     enddo
   enddo
  end subroutine
!----------------------------------------------------------------------------------------------------------------------------
!  subroutine contact_nrg(xyz,nrg,cmap,nres)
!   !---------------------------------------
!   implicit none
!   !---------------------------------------
!   !real,parameter::contact_cutsq=contact_cut**2
!   !---------------------------------------
!   integer,intent(in)::nres
!   real,intent(out)::nrg
!   real,dimension(:,:),intent(inout)::xyz
!   real,dimension(:,:),intent(in)::cmap
!   !---------------------------------------
!   integer::ires
!   integer::jres
!   real::x
!   real::dij
!   real::dijsq
!   real,dimension(3)::rij
!   !---------------------------------------
!   nrg=0.
!   do ires = 1, nres-nlocal
!     do jres = ires + nlocal, nres
!       if(cmap(ires,jres) /= 0) then
!         rij = xyz(:,ires) - xyz(:,jres)
!         dijsq = sum(rij*rij)
!         !if(dijsq <= contact_cutsq) then
!           dij = sqrt(dijsq)
!           x   = k_contact*(0.7*(atan(x-8.0))/pi-0.5)
!           !x   = k_contact*((atan(x-8.0))/pi-0.5)
!           nrg = nrg + x
!         !endif
!       endif
!     enddo
!   enddo
!   !---------------------------------------
!  end subroutine contact_nrg
!!----------------------------------------------------------------------------------------------------------------------------
!  subroutine contact_force(xyz,cmap,contactf,nres)
!   !---------------------------------------
!   implicit none
!   !---------------------------------------
!   real,parameter::contact_cutsq=contact_cut**2
!   !---------------------------------------
!   integer,intent(in)::nres
!   real,dimension(:,:),intent(in)::cmap
!   real,dimension(:,:),intent(inout)::xyz
!   real,dimension(:,:),intent(inout)::contactf
!   !---------------------------------------
!   integer::ires
!   integer::jres
!   real::x
!   real::dij
!   real::dijsq
!   real,dimension(3)::rij
!   !---------------------------------------
!   contactf=0
!   do ires = 1, nres-nlocal
!     do jres = ires + nlocal, nres
!       if(cmap(ires,jres) /= 0) then
!         rij   = xyz(:,ires) - xyz(:,jres)
!         dijsq = sum(rij*rij)
!         if(dijsq <= contact_cutsq) then
!           dij   = sqrt(dijsq)
!           x     = -k_contact*(0.7/pi/(0.7*(x-8.0)**2+1.0))
!           !x     = -k_contact*(1.0/pi/((x-8.0)**2+1.0))
!           contactf(:,ires) = contactf(:,ires) + x*rij/dij    
!           contactf(:,jres) = contactf(:,jres) - x*rij/dij 
!         endif
!       endif
!     enddo
!   enddo
!   !---------------------------------------
!  end subroutine contact_force
!----------------------------------------------------------------------------------------------------------------------------
  subroutine berendsen_tcoup(t,curr_t,vel,nres)
   !---------------------------------------
   implicit none
   !---------------------------------------
   integer,intent(in)::nres
   real,intent(in)::t
   real,dimension(:,:),intent(inout)::vel
   !---------------------------------------
   integer::ires
   real::ke
   real::lamba
   real::curr_t !I am assuming intent(in)
   !---------------------------------------
   ke = 0.0
   do ires=1,nres
      ke = ke + sum(vel(:,ires)*vel(:,ires)) 
   enddo
   ke = ke*0.5
   curr_temp = ke/real((3*nres-3)/2.)
   lamba = sqrt(1+tcoup*(t/curr_t-1))
   vel = vel*lamba
   !---------------------------------------
  end subroutine 
!----------------------------------------------------------------------------------------------------------------------------
  subroutine random_force(model,ran_f,t,zi,dt,nres)
   !---------------------------------------
   implicit none
   !---------------------------------------
   integer,intent(in)::nres
   real,intent(in)::t
   real,intent(in)::zi
   real,intent(in)::dt
   real,dimension(:,:),pointer::model
   real,dimension(:,:),pointer::ran_f
   !---------------------------------------
   integer::jres
   integer::i
   real::x
   !---------------------------------------
   ran_f=0.
   do jres=1,nres
      do i=1,3
         if (jres/=go_s .and. jres/=go_e) then !to define region of go app
         ran_f(i,jres)=ran_gauss()
         endif
      enddo
   enddo
   x=sqrt(2.0*t*zi/dt)
   ran_f=ran_f*x
   !---------------------------------------
  end subroutine
!----------------------------------------------------------------------------------------------------------------------------
  subroutine dihedral_nrg(xyz,nrg,nres,targt)
   !---------------------------------------
   implicit none
   !---------------------------------------
   integer,intent(in)::nres
   real,intent(out)::nrg
   real,dimension(:),intent(in)::targt
   real,dimension(:,:),intent(in)::xyz
   !---------------------------------------
   integer::ires
   real::x
   real::ao
   real::a
   real,dimension(nres)::angles
   !---------------------------------------
   nrg=0.
   call update_dihedral(angles,xyz,nres)
   do ires=2,nres-2
     a =angles(ires)*deg2rad
     ao=targt(ires)*deg2rad
     !x = k_dih1*(1+cos(a-ao))+ k_dih3*(1+cos(3*(a-ao)))
     x = (k_dih1*(1+cos(a-ao)) + k_dih3*(1+cos(3*(a-ao))))**expo/((k_dih1+k_dih3)*2)**(expo-1)
     nrg = nrg + x
   enddo
   !---------------------------------------
  end subroutine 
!----------------------------------------------------------------------------------------------------------------------------
  subroutine dihedral_force(xyz,dih_f,nres,targt)
   !---------------------------------------
   implicit none
   !---------------------------------------
   integer,intent(in)::nres
   real,dimension(:,:),intent(in)::xyz
   real,dimension(:,:),pointer::dih_f
   real,dimension(:),pointer::targt
   !---------------------------------------
   integer::ires
   real::x
   real::v1,v2,v3
   real::a
   real::ao
   real::angle
   real,dimension(3,4)::frame
   real,dimension(3,4)::gradient
   !---------------------------------------
   dih_f=0.
   do ires=2,nres-2
     frame = xyz(:,ires-1:ires+2)
     call dihedral_gradient(frame,gradient,angle)
     a =angle*deg2rad
     ao=targt(ires)*deg2rad
     !x = k_dih*(sin(a-ao)) + k_dih3*3*(sin(3*(a-ao)))
     !x = (expo*(k_dih*(1+cos(a-ao))+k_dih3*(1+cos(3*(a-ao))))**(expo-1))*(k_dih*sin(a-ao) + k_dih3*3*sin(3*(a-ao))) 
     v1 = (expo*(k_dih1*(1+cos(a-ao))+k_dih3*(1+cos(3*(a-ao)))))
     v2 = v1**(expo-1)
     v3 = (k_dih1*sin(a-ao) + k_dih3*3*sin(3*(a-ao)))
     !write(*,*) v1,v2,v3 
     !/((k_dih+k_dih3)*2)**(expo-1)
     x = v3*v2
     !if(ires/=go_s .and. ires/=go_e) then
     dih_f(:,ires-1) = dih_f(:,ires-1) + x*gradient(:,1)
     dih_f(:,ires)   = dih_f(:,ires)   + x*gradient(:,2)
     dih_f(:,ires+1) = dih_f(:,ires+1) + x*gradient(:,3)
     dih_f(:,ires+2) = dih_f(:,ires+2) + x*gradient(:,4)
     !endif
   enddo
   !---------------------------------------
  end subroutine 
!----------------------------------------------------------------------------------------------------------------------------
  subroutine anglebend_nrg(xyz,nrg,nres,targt)
   !---------------------------------------
   implicit none
   !---------------------------------------
   integer,intent(in)::nres
   real,intent(out)::nrg
   real,dimension(:),pointer::targt
   real,dimension(:,:),intent(in)::xyz
   !---------------------------------------
   integer::ires
   real::x
   real::a
   real::ao
   real,dimension(nres)::angles
   !---------------------------------------
   nrg=0.
   call update_anglebend(angles,xyz,nres)
   do ires=2,nres-1
     a =angles(ires)*deg2rad
     ao=targt(ires)*deg2rad
     x = k_anglebend*0.5*(a-ao)**2
     nrg =nrg + x
   enddo
   !---------------------------------------
  end subroutine 
!----------------------------------------------------------------------------------------------------------------------------
  subroutine anglebend_force(xyz,bend_f,nres,targt)
   !---------------------------------------
   implicit none
   !---------------------------------------
   integer,intent(in)::nres
   real,dimension(:),pointer::targt
   real,dimension(:,:),intent(in)::xyz
   real,dimension(:,:),pointer::bend_f
   !---------------------------------------
   integer::ires
   real::angle
   real::x
   real::a
   real::ao
   real,dimension(3,3)::frame
   real,dimension(3,3)::gradient
   !---------------------------------------
   bend_f=0.
   do ires=2,nres-1
     frame = xyz(:,ires-1:ires+1)
     call anglebend_gradient(frame,gradient,angle)
     a =angle*deg2rad
     ao=targt(ires)*deg2rad
     x = - k_anglebend*(a-ao)
     !if(ires/=go_s .and. ires/=go_e) then
     bend_f(:,ires-1) = bend_f(:,ires-1) + x*gradient(:,1) 
     bend_f(:,ires+1) = bend_f(:,ires+1) + x*gradient(:,3) 
     bend_f(:,ires)   = bend_f(:,ires)   + x*gradient(:,2)
     !endif 
   enddo
   !---------------------------------------
  end subroutine 
!----------------------------------------------------------------------------------------------------------------------------
  subroutine update_dihedral(angs,xyz,nres)
   !---------------------------------------
   implicit none
   !---------------------------------------
   integer,intent(in)::nres
   real,dimension(:,:),intent(in)::xyz
   real,dimension(:),intent(inout)::angs
   !---------------------------------------
   integer::ires
   real::sb
   real::sb2
   real::sa_c
   real::sa_b
   real::sb_c
   real::d
   real::e
   real::g
   real::theta
   real,dimension(3)::va
   real,dimension(3)::vb
   real,dimension(3)::vc
   real,dimension(3)::c_vb_vc
   real,dimension(3)::ub
   real,dimension(3)::iva
   !---------------------------------------
   angs=0
   do ires=1,nres-3 !why it is done like this is beyond me 
     va  = xyz(:,ires+1)-xyz(:,ires)   
     vb  = xyz(:,ires+2)-xyz(:,ires+1)  
     vc  = xyz(:,ires+3)-xyz(:,ires+2) 
     sb2 = sum(vb*vb)
     sb  = sqrt(sb2)
     call crossprod(vb,vc,c_vb_vc) 
     sa_c = sum(-va*vc)
     sa_b = sum(va*vb)
     sb_c = sum(vb*vc)
     e =  sum(va*c_vb_vc)/sb
     g =  sa_c + sa_b*sb_c/sb2
     theta = atan2(e,g)*rad2deg
     if(theta<0.) theta=theta+360.
     angs(ires+1)=theta
     !write(*,*)ires+1,angs(ires+1)
   enddo
   !---------------------------------------
  end subroutine update_dihedral
   !---------------------------------------
   !----------------------------------------------------------------------------------------------------------------------------
  subroutine update_anglebend(angs,xyz,nres)
   !---------------------------------------
   implicit none
   !---------------------------------------
   integer,intent(in)::nres
   real,dimension(:,:),intent(in)::xyz
   real,dimension(nres),intent(inout)::angs
   !---------------------------------------
   integer::ires
   real::angle
   real::d21
   real::d23
   real,dimension(3)::v21
   real,dimension(3)::v23
   !---------------------------------------
    angs = 0.
    do ires=2,nres-1     
      v21 = xyz(1:3,ires-1) - xyz(1:3,ires)
      v23 = xyz(1:3,ires+1) - xyz(1:3,ires)
      d21 = sqrt(sum(v21*v21))
      d23 = sqrt(sum(v23*v23))
      angle = (sum(v21*v23))/(d21*d23)
      if(angle.gt.1)  angle = 1
      if(angle.lt.-1) angle = -1
      angle = acos(angle)*rad2deg
      angs(ires) = angle
      !write(*,*)angs(ires)
    enddo
   !---------------------------------------
  end subroutine update_anglebend
!----------------------------------------------------------------------------------------------------------------------------
  subroutine get_cmap(protein,nres,cmap) 
   !------------------------
   implicit none
   !------------------------
   real,parameter::allatm_tolsq=allatm_tol*allatm_tol
   !real,parameter::calpha_tolsq=calpha_tol*calpha_tol
   !------------------------
   integer,intent(inout)::nres
   type(residue),dimension(:),pointer::protein
   real,dimension(:,:),pointer::cmap
   !------------------------
   integer::c
   integer::iatom
   integer::jatom
   integer::ires
   integer::jres
   integer::ncon
   integer::ncon2
   real::rijsq
   real::rij
   !------------------------------------------------------------------------------------
   ncon=0
   allocate(cmap(nres,nres))
   cmap=0
   !do ires=1,nres-nlocal
   !  do jres=ires+nlocal,nres
   !     rijsq=sum((protein(ires)%xyz(:,2)-protein(jres)%xyz(:,2))**2)
   !     if(rijsq <= calpha_tolsq) then
   !       rij = sqrt(rijsq)
   !       cmap(ires,jres) = rij
   !       cmap(jres,ires) = rij**10
   !       ncon=ncon+1
   !     endif
   !  enddo
   !enddo
   !------------------------------------------------------------------------------------
   !write(*,*)ncon
   !ncon=0
   !ncon2=0
   !------------------------------------------------------------------------------------
   do ires=1,nres-nlocal
     do jres=ires+nlocal,nres
       loop3:do iatom=1,22
         do jatom=1,22
           if(protein(ires)%xyz(1,iatom) /= 999.0 .and. protein(jres)%xyz(1,jatom) /= 999.0) then 
             rijsq=sum((protein(ires)%xyz(:,iatom) - protein(jres)%xyz(:,jatom))**2)
             if(rijsq <= allatm_tolsq) then
               rijsq=sum((protein(ires)%xyz(:,2)-protein(jres)%xyz(:,2))**2)
!               if(rijsq <= calpha_tolsq) then
                 rij = sqrt(rijsq)
                 cmap(jres,ires) = rij**10
                 cmap(ires,jres) = rij
                 ncon=ncon+1
                 exit loop3
!               !else
!               !  write(*,*)ires,jres,iatom,jatom,sqrt(sum((protein(ires)%xyz(:,2)-protein(jres)%xyz(:,2))**2))
!               endif
!               !ncon2=ncon2+1
             endif
           endif
        enddo
       enddo loop3
     enddo
   enddo
   !write(*,*)ncon,ncon2
   !------------------------------------------------------------------------------------
!   do ires=0,nres
!      write(*,'(i2,$)')ires
!   enddo
!   write(*,*)
!   do ires=1,nres
!     write(*,'(i2,$)')ires
!     do jres=1,nres
!       c=0
!       if(cmap(ires,jres)/=0) c=1
!       write(*,'(i2,$)')c
!     enddo
!     write(*,*)
!   enddo
   !------------------------------------------------------------------------------------
  end subroutine
!----------------------------------------------------------------------------------------------------------------------------
  subroutine init_ranseed(process_id)
   !------------------------------------------
   implicit none
   !------------------------------------------
   integer,dimension(8) :: x
   integer,intent(in)::process_id
   integer::isize,i
   integer,dimension(:),allocatable::iseed
   !------------------------------------------
   call date_and_time(values=x)
   call random_seed(size=isize)
   allocate(iseed(isize))
   call random_seed(get=iseed)
   iseed = (iseed +process_id) * (x(8)-500)
   call random_seed(put=iseed)
   deallocate(iseed)
   !------------------------------------------
  end subroutine
!----------------------------------------------------------------------------------------------------------------------------
  subroutine write_traj(iunit,xyz,nres,istep,seq,temp)
   !------------------------------------------
   implicit none
   !------------------------------------------
   character(len=3),dimension(20),parameter::three=(/'ALA','CYS','ASP','GLU','PHE','GLY','HIS','ILE',&
     'LYS','LEU','MET','ASN','PRO','GLN','ARG','SER','THR','VAL','TRP','TYR'/)
   !------------------------------------------
   integer,intent(in)::nres
   integer,intent(in)::iunit
   integer,intent(in)::istep
   integer,dimension(:),intent(in)::seq
   real,intent(in)::temp
   real,dimension(:,:),pointer::xyz
   !------------------------------------------
   integer::ires
   integer::irec=0
   !------------------------------------------
   if(iunit/=99) then
     do ires=1,nres
       write(iunit,'("ATOM",3X,i4,2X,"CA  ",a3,2x,i4,4X,3F8.3)')ires,three(iseq(ires)),ires,xyz(1:3,ires) 
     enddo
     write(iunit,'("TER")')
     write(iunit,'("ENDMDL",i12)')istep
   else
     irec=irec+1
     write(iunit,rec=irec)temp,istep,xyz(1:3,1:nres)
   endif
   !do ires=1,nres
   !   write(iunit,'("ATOM",3X,i4,2X,"CA  ",a3,2X,i4,4X,3F8.3)')ires,three(seq(ires)),ires,xyz(1:3,ires)
   ! enddo
   ! write(iunit,'("TER",i10)')step
   ! write(iunit,'("ENDMDL")')
   !------------------------------------------
  end subroutine
!----------------------------------------------------------------------------------------------------------------------------
  real function ran_gauss()
   !------------------------------------
   implicit none
   !------------------------------------
   integer::iset=0   
   real::fac
   real::gset
   real::rsq
   real::v1
   real::v2
   real::ran1
   real::x
   real::y
   save::iset  
   save::gset  
   !------------------------------------
   rsq=2.0
   if(iset.eq.0) then
    do while(rsq >= 1.0 .or. rsq <= 0.)
      call random_number(x)
      call random_number(y)
      v1=2.0*x-1
      v2=2.0*y-1
      rsq=v1*v1+v2*v2
    enddo
      fac=sqrt(-2.0*log(rsq)/rsq)
      gset=v1*fac
      ran_gauss=v2*fac
      iset=1
   else
      ran_gauss=gset
      iset=0
   endif
   !------------------------------------
  end function ran_gauss
!----------------------------------------------------------------------------------------------------------------------------
  subroutine bond_nrg(xyz,nres,nrg)
   !------------------------------------------
   implicit none
   !------------------------------------------
   integer,intent(in)::nres
   real,intent(inout)::nrg
   real,dimension(:,:),intent(inout)::xyz
   !------------------------------------------
   integer::ires
   integer::jres
   real::dij
   real::x
   real,dimension(3)::rij
   !------------------------------------------
   nrg=0.
   do ires = 1, nres-1
     jres = ires + 1 
     rij = xyz(:,ires) - xyz(:,jres)
     dij = sqrt(sum(rij*rij))
     nrg =  nrg + k_bond*0.5*(dij-3.8)**2
   enddo
   !------------------------------------------
  end subroutine
!-----------------------------------------------------------------------------------------------------------------
  subroutine dsbond_nrg(dsbonds,native_s,xyz,nres,ndsbonds,nrg)
   !------------------------------------------
   implicit none
   !------------------------------------------
   integer,dimension(:,:),pointer::dsbonds
   integer::ndsbonds
   real,dimension(:,:)::native_s
   integer,intent(in)::nres
   real,intent(inout)::nrg
   real,dimension(:,:),intent(in)::xyz
   !------------------------------------------
   integer::ires
   integer::idsbonds
   integer::jres
   real::d ! current distance between the two alpha carbons
   real::doo ! native distance between the two alpha carbons
   real::x
   real,dimension(3)::r ! displacement vector between the two alpha carbons
   real,dimension(3)::ro ! native displacement vector between the two alpha carbons 
   !------------------------------------------
   nrg=0.
   do ires = 1, nres 
     do idsbonds=1,ndsbonds
       if(ires == dsbonds(1,idsbonds)) then
          r = xyz(:,dsbonds(1,idsbonds)) - xyz(:,dsbonds(2,idsbonds))
          ro = native_s(:,dsbonds(1,idsbonds)) - native_s(:,dsbonds(2,idsbonds))
          d = sqrt(sum(r*r))
          doo = sqrt(sum(ro*ro))
          nrg =  nrg + k_bond*0.5*(d-doo)**2
       endif
     end do
   enddo
   !------------------------------------------
  end subroutine
!----------------------------------------------------------------------------------------------------------------------------
   subroutine sdhmon(ange,xyz,dsbonds) !disulfide dihedral monitor
   !---------------------------------------
   implicit none
   !---------------------------------------
   integer,dimension(:,:),pointer::dsbonds
   real,dimension(:,:),intent(in)::xyz
   real,dimension(:),pointer::ange
   !---------------------------------------
   integer::ires
   integer::jres
   real::sb
   real::sb2
   real::sa_c
   real::sa_b
   real::sb_c
   !real::d
   real::e
   real::g
   real::theta
   real,dimension(3)::va
   real,dimension(3)::vb
   real,dimension(3)::vc
   real,dimension(3)::c_vb_vc
   !real,dimension(3)::ub
   !real,dimension(3)::iva
   !---------------------------------------
   ange=0
   do ires=1,ndsbonds
     va  = xyz(:,dsbonds(1,ires)+1)-xyz(:,dsbonds(1,ires))   
     vb  = xyz(:,dsbonds(1,ires))-xyz(:,dsbonds(2,ires))  
     vc  = xyz(:,dsbonds(2,ires))-xyz(:,dsbonds(2,ires)-1) 
     sb2 = sum(vb*vb)
     sb  = sqrt(sb2)
     call crossprod(vb,vc,c_vb_vc) 
     sa_c = sum(-va*vc)
     sa_b = sum(va*vb)
     sb_c = sum(vb*vc)
     e =  sum(va*c_vb_vc)/sb
     g =  sa_c + sa_b*sb_c/sb2
     theta = atan2(e,g)*rad2deg
     !if(theta<0.) theta=theta+360.
     if(theta<-180) then
       theta=theta+360
     else if(theta>180)then
       theta=theta-360
     endif
     ange(ires)=theta
   enddo
   !do jres=1,ndsbonds
     !write(*,*),"Step: ",istep," Disulfide dihedral",dsbonds(1,jres), &
    !     "-",dsbonds(2,jres),ange(jres),";"
   !enddo
   !deallocate(ange)
   !---------------------------------------
  end subroutine sdhmon
!----------------------------------------------------------------------------------------------------------------------------
  subroutine sdhdir(angpr,angdr,rotdr,sumrr)
   !------------------------------------------
   implicit none
   !------------------------------------------
   real,dimension(:),pointer::angpr
   real,dimension(:),pointer::angdr
   real,dimension(:),pointer::rotdr
   real,dimension(:),pointer::sumrr
   !-----------------------------------------
   integer::i !loop counter
   real :: prev,curr
   !write(*,*) 'Rotd calculated'
   !rotdr=angdr-angpr
   !write(*,*) 'ROTD:',rotdr(:)
   !write(*,*)'NDSBONDS-----',ndsbonds
   do i=1,ndsbonds
     if(.not. associated(rotdr)) then
        write(*,*)'Lost rotd'
        exit
     endif
     !write(*,*)'rotd',rotdr(:) 
     prev=angpr(i)
     curr=angdr(i)
     !if(prev<0) prev=prev+360 !convert negatives 
     !if(curr<0) curr=curr+360 
     rotdr(i)=curr-prev !difference between current and previous ds dihedral
     if (rotdr(i)>180 ) then
        rotdr(i)=360-rotdr(i)
     else if(rotdr(i)<-180) then
        rotdr(i)=-360-rotdr(i)
     endif
     if (rotdr(i) >0) then
        !write(*,*)'before'
        sumrr(i)=sumrr(i)+abs(rotdr(i)) !anti-clokwise flag wrt corkscrew rule
        !write(*,*)'Rotd ',rotdr(i)
     else if(rotdr(i)<0) then
         !write(*,*)'before'
        sumrr(i)=sumrr(i)-abs(rotdr(i)) !clockwise flag wrt corkscrew rule
        !write(*,*)'Rotd',rotdr(:)
     endif
   enddo
   angpr=angdr 
   !write (*,*)'sumr test',sumrr(:)
   end subroutine sdhdir
!----------------------------------------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------------------------------------  
  subroutine bond_force(xyz,nres,bondf)
   !------------------------------------------
   implicit none
   !------------------------------------------
   integer,intent(in)::nres
   real,dimension(:,:),intent(inout)::xyz
   real,dimension(:,:),pointer::bondf
   !------------------------------------------
   integer::ires
   integer::jres
   real::dij
   real::x
   real,dimension(3)::rij
   !------------------------------------------
   bondf=0.
   do ires = 1, nres-1
     jres = ires + 1 
     rij = xyz(:,ires) - xyz(:,jres)
     dij = sqrt(sum(rij*rij))
     x =  - k_bond*(dij-3.8)
     !if(ires/=go_s .and. ires/=go_e) then
     bondf(:,ires) = bondf(:,ires) + x*rij/dij    
     bondf(:,jres) = bondf(:,jres) - x*(rij/dij)
     !endif 
   enddo
   !------------------------------------------
  end subroutine bond_force
  subroutine dsbond_force(xyz,native_s,dsbonds,ndsbonds,nres,dsbondf)
   !------------------------------------------
   implicit none
   !------------------------------------------
   integer,intent(in)::nres
   integer,intent(in)::ndsbonds
   integer,pointer,dimension(:,:)::dsbonds
   real,dimension(:,:)::native_s
   real,dimension(:,:),intent(inout)::xyz
   real,dimension(:,:),pointer::dsbondf
   !------------------------------------------
   integer::ires
   integer::jres
   logical::isInDSBond
   real::d
   real::doo
   real::x
   integer::idsbonds
   real,dimension(3)::r
   real,dimension(3)::ro
   !------------------------------------------
   dsbondf=0.
   isInDSBond = .false.
  ! do ires = 1, nres
     do idsbonds=1,ndsbonds
     !  if(ires == dsbonds(1,idsbonds)) then !if the current residue is found in the disulfide bond matrix
         r = xyz(:,dsbonds(1,idsbonds)) - xyz(:,dsbonds(2,idsbonds)) ! calculate current disulfide bond displacement vector
         ro = native_s(:,dsbonds(1,idsbonds)) - native_s(:,dsbonds(2,idsbonds)) ! calculate native disulfide bond displacement vector
         d = sqrt(sum(r*r)) !find r vector magnitiude
         doo= sqrt(sum(ro*ro)) !find ro vector maginitude
         x =  - k_bond*(d-doo) !calculate magnitude of the bond force
         !dsbondf(:,dsbonds(1,idsbonds)) = dsbondf(:,dsbonds(idsbonds,1)) + x*r/d    
         !dsbondf(:,dsbonds(2,idsbonds)) = dsbondf(:,dsbonds(idsbonds,2)) - x*(r/d)
         dsbondf(:,dsbonds(1,idsbonds)) = +x*r/d
         dsbondf(:,dsbonds(2,idsbonds)) = -x*(r/d)

     !    endif
     enddo 
   !enddo
   !------------------------------------------
  end subroutine dsbond_force
!----------------------------------------------------------------------------------------------------------------------------
  subroutine rand_chain(xyz,nres)
   !------------------------------------------
   implicit none
   !------------------------------------------
   integer,intent(in)::nres
   real,dimension(3,nres),intent(inout)::xyz
   !------------------------------------------
   integer::ires
   integer::ios
   real::randx
   real::x,y,z,angle,delta,xy
   real,dimension(3)::v21,v23
   real::d21,d23
   logical::clashing
   !------------------------------------------
   angle=0; xyz=0
   xyz(1:3,1) = (/0.,0.,0./)
   xyz(1:3,2) = (/3.8,0.,0./)
   do while(angle*rad2deg < 70 .or. angle*rad2deg > 140)
     call random_number(randx)
     angle = randx*pi   
   enddo
   y = sin(angle)*3.8
   x = 3.8 + cos(angle)*3.8
   xyz(1:3,3) = (/x,y,0.0/)
   !!-----------------------------------------
   ires=4
   do 
     clashing=.false.
     call random_number(randx)
     angle = randx*pi 
     z = sin(angle)*3.8
     if(angle > (pi/2.)) then 
       angle = - angle + pi/2.
       z = sin(angle)*3.8
     endif
     z = xyz(3,ires-1) + z
     xy = cos(angle)*3.8
     call random_number(randx)
     angle = (randx*pi)  
     y = sin(angle)*xy
     x = cos(angle)*xy
     y = xyz(2,ires-1) + y
     x = xyz(1,ires-1) + x
     xyz(1:3,ires)=(/x,y,z/)
     !---------- Check Delta ---------------
     v21 = xyz(1:3,ires-2) - xyz(1:3,ires-1)
     v23 = xyz(1:3,ires) - xyz(1:3,ires-1)
     d21 = sqrt(sum(v21*v21))
     d23 = sqrt(sum(v23*v23))
     angle = (sum(v21*v23))/(d21*d23)
     if(angle > 1)  angle = 1; if(angle < -1) angle = -1
     angle = acos(angle)*rad2deg
     delta = angle
     !--------------------------------------
     clashing=clash(xyz,ires,nres)
     if(delta > 70 .and. delta < 140 .and. .not.clashing) ires = ires + 1 
     if(ires==nres+1) exit
   enddo
   !do ires=1,nres
   !  write(*,'("ATOM",3X,i4,2X,"CA  ","ALA",2X,i4,4X,3F8.3)')ires,ires,xyz(1:3,ires)
   !enddo
   !------------------------------------------
  end subroutine rand_chain
!----------------------------------------------------------------------------------------------------------------------------
  logical function clash(xyz,ires,nres)
   !------------------------------------------
   implicit none
   !------------------------------------------
   integer,intent(in)::ires
   integer,intent(in)::nres
   real,dimension(:,:),intent(in)::xyz
   !------------------------------------------
   integer::jres
   real::d
   real::x
   real::y
   real,dimension(3)::v
   !------------------------------------------
   clash=.false.
   do jres=1,nres
     if(jres==ires)cycle
     v = xyz(:,jres) - xyz(:,ires)
     d = sqrt(sum(v*v))
     if(d < 3.5) then
      clash=.true.
     endif
   enddo
   !------------------------------------------
  end function clash
!----------------------------------------------------------------------------------------------------------------------------
  subroutine get_rmsd(c0,c1,err,natm)
   !------------------------------------------
   ! coord0 is template. coord1 is target. err is rmsd of natm atoms
   !------------------------------------------
   implicit none
   !------------------------------------------
   real,parameter::tol=0.00001
   !------------------------------------------
   integer,intent(in)::natm
   real,intent(out)::err
   real,dimension(3,natm),intent(in)::c0
   real,dimension(3,natm),intent(in)::c1
   !------------------------------------------
   integer::ix
   integer::iy
   integer::iz
   integer::i
   integer::j
   integer::k
   integer::iflag
   integer::ict
   integer::nnn
   integer::ios
   real::y
   real,dimension(3)::xc0
   real,dimension(3)::xc1
   real(8)::sig
   real(8)::sg
   real(8)::gam
   real(8)::bet
   real(8)::cc
   real(8)::bb
   real(8),dimension(3)::t
   real(8),dimension(3,3)::aa
   real(8),dimension(3,3)::rot
   real,dimension(:,:),allocatable::x0,x1
   !------------------------------------------
   if (allocated(x0)) deallocate(x0)
   allocate(x0(3,natm),stat=ios)
   if (ios/=0) stop 'subroutine get_rmsd :: allocation failed for x0'
   if (allocated(x1)) deallocate(x1)
   allocate(x1(3,natm),stat=ios)
   if (ios/=0) stop 'subroutine get_rmsd ::  allocation failed for x1'
   xc0 = 0.
   xc1 = 0.
   do I=1,natm
     xc0 = xc0 + c0(1:3,I)
     xc1 = xc1 + c1(1:3,I)
   enddo
   if (natm < 4) then
     write(*,'("subroutine get_rmsd ::  Number of atoms   =",i9,"  <===too small.")') natm
     err = -1
     return
   endif
   xc0 = xc0/real(natm)   !! center of mass, target
   xc1 = xc1/real(natm)   !! center of mass, template
   do I=1,natm
     x0(:,I) = c0(:,I) - xc0
     x1(:,I) = c1(:,I) - xc1
   enddo
   !* x0 and x1 are centered on the origin
   aa = 0.
   do i=1,3
     do j=1,3
       aa(i,j) = dot_product(x1(i,:),x0(j,:))
     enddo
   enddo
   rot = 0.
   do i=1,3
     rot(i,i) = 1.
   enddo
   ict=0
   iflag=0
   ix=1
   do while (ict < 1000)
     ict=ict+1
     do ix=1,3
       iy=mod(ix,3)+1
       iz=6-ix-iy
       sig=aa(iz,iy)-aa(iy,iz)
       gam=aa(iy,iy)+aa(iz,iz)
       sg=sqrt(sig*sig+gam*gam)
       if(sg == 0.) cycle
       sg=1./sg
       if(abs(sig) <= tol*abs(gam)) cycle
       do k=1,3
         bb=gam*aa(iy,k)+sig*aa(iz,k)
         cc=gam*aa(iz,k)-sig*aa(iy,k)
         aa(iy,k)=bb*sg
         aa(iz,k)=cc*sg
         bb=gam*rot(iy,k)+sig*rot(iz,k)
         cc=gam*rot(iz,k)-sig*rot(iy,k)
         rot(iy,k)=bb*sg
         rot(iz,k)=cc*sg
       enddo
       iflag=1
     enddo
     if (iflag==0) exit
   enddo
   err=0.
   do j=1,3
     t(j) = dot_product(rot(j,:),xc1(:))
   enddo
   do I=1,natm
     do j=1,3
       t(j)=dot_product(rot(j,:),x1(:,I))
     enddo
     x1(:,I)= t + xc0  !! x1 is now superimposed on x0
     t = c0(:,I) - x1(:,I)
     err = err + dot_product(t,t)
   enddo
   err=sqrt(err/natm)
   deallocate(x0,x1)
   !------------------------------------------
  end subroutine get_rmsd
!----------------------------------------------------------------------------------------------------------------------------
  subroutine get_Q(xyz,cmap,nres,q) 
   !------------------------------------------
   implicit none
   !------------------------------------------
   integer,intent(in)::nres
   real,intent(out)::q
   real,dimension(:,:),intent(in)::xyz
   real,dimension(:,:),pointer::cmap
   !------------------------------------------
   integer::jres,ncont
   real::dsq
   real,dimension(3)::rij
   !------------------------------------------
   q=0; ncont=0
   do ires=1,nres
     do jres=ires+3,nres
       if(cmap(ires,jres)/=0)then
          ncont=ncont+1
          rij = xyz(:,ires) - xyz(:,jres)
          dsq = sqrt(sum(rij*rij))
          if(dsq < cmap(ires,jres)*1.2)then
             q = q + 1
          endif
       endif
     enddo
   enddo
   if(ncont/=0) then
     q=q/real(ncont)
   else
     q=999.
   endif
   !------------------------------------------
  end subroutine get_Q
!----------------------------------------------------------------------------------------------------------------------------
  subroutine anglebend_gradient(xyz,grad,delta)
   !--------------------------------------
   implicit none
   !--------------------------------------
   real,intent(out)::delta
   real,dimension(3,3),intent(in)::xyz
   real,dimension(3,3),intent(out)::grad
   !--------------------------------------
   integer::i,j
   real::a
   real::b
   real::a2
   real::b2
   real::ab
   real::angle
   real::isinang
   real::cosang
   real,dimension(3)::va
   real,dimension(3)::vb
   real,dimension(3)::vab
   !--------------------------------------
   grad=0
   va = xyz(1:3,1) - xyz(1:3,2)
   vb = xyz(1:3,3) - xyz(1:3,2)
   vab = va*vb
   a2 = sum(va*va)
   a  = sqrt(a2)
   b2 = sum(vb*vb)
   b  = sqrt(b2)
   ab = a*b
   angle = sum(vab)/(a*b)
   if(angle>1) angle = 1; if(angle<-1) angle = -1
   angle = acos(angle)
   isinang= -1/sin(angle)
   cosang=cos(angle)
   do i=1,3 !! X,Y,Z     K(1) <--(a)-- L(2) ---(b)---> M(3)
     grad(i,1)= isinang*(cosang * ( xyz(i,2)-xyz(i,1) )/a2 + ( xyz(i,3)-xyz(i,2) )/ab )
     grad(i,3)= isinang*(cosang * ( xyz(i,2)-xyz(i,3) )/b2 + ( xyz(i,1)-xyz(i,2) )/ab )
   enddo
   grad(:,2)= - grad(:,1) - grad(:,3)
   delta=angle*rad2deg
   !--------------------------------------
  end subroutine anglebend_gradient
!----------------------------------------------------------------------------------------------------------------------------
  subroutine dihedral_gradient(xyz,grad,theta)
   !--------------------------------------
   implicit none
   !--------------------------------------
   real,intent(out)::theta
   real,dimension(3,4),intent(in)::xyz
   real,dimension(3,4),intent(out)::grad
   !--------------------------------------
   integer::i,jitr
   real::sb
   real::sb2
   real::sa_c
   real::sa_b
   real::sb_c
   real::d
   real::e
   real::f
   real::g
   real::h
   real::j
   real::k
   real::l
   real::m
   real::dtheta_dd
   real,dimension(3)::va
   real,dimension(3)::vc
   real,dimension(3)::c_vb_vc
   real,dimension(3)::ub
   real,dimension(3)::iva
   real,dimension(3)::vb
   real,dimension(4)::dh_dx
   real,dimension(4)::dk_dx
   real,dimension(4)::dj_dx
   real,dimension(4)::dl_dx
   real,dimension(4)::de_dx
   real,dimension(4)::dm_dx
   real,dimension(4)::dd_dx
   logical :: nan
   !--------------------------------------
   grad=0
   nan = .false.
   !--------------------------------------
   va  = xyz(:,2)-xyz(:,1) !v21
   vb  = xyz(:,3)-xyz(:,2) !v32 
   vc  = xyz(:,4)-xyz(:,3) !v43
   sb2 = sum(vb*vb)
   sb  = sqrt(sb2)
   call crossprod(vb,vc,c_vb_vc)
   sa_c = sum(-va*vc)
   sa_b = sum(va*vb)
   sb_c = sum(vb*vc)
   e =  sum(va*c_vb_vc)/sb
   g =  sa_c + sa_b*sb_c/sb2
   f = 1./g
   h = sum(-va*vc)
   j = sum(va*vb)
   k = sum(vb*vc)
   l = 1/sb2
   m = e*sb 
   if(g == 0)return
   d = e/g
   theta = atan2(e,g)*rad2deg
   if(theta<0) theta=theta+360.
   dtheta_dd = 1/(1.+d*d)
   do i=1,3   !! X,Y,Z   1 ----(a)---> 2 ----(b)----> 3 ----(c)---> 4
       dh_dx(1)= vc(i) ; dh_dx(2)= -vc(i)        ; dh_dx(3)= va(i)      ; dh_dx(4)= -va(i);
       dk_dx(1)= 0     ; dk_dx(2)= -vc(i)        ; dk_dx(3)= vc(i)-vb(i); dk_dx(4)= vb(i);
       dj_dx(1)= -vb(i); dj_dx(2)=  vb(i)-va(i)  ; dj_dx(3)= va(i)      ; dj_dx(4)= 0;
       dl_dx(2)=  2*(xyz(i,3)-xyz(i,2))/(sb2*sb2); dl_dx(1)= 0 
       dl_dx(3)= -2*(xyz(i,3)-xyz(i,2))/(sb2*sb2); dl_dx(4)= 0  
       !!----------------------------------------------------------------------------
       if(i==1)then
         dm_dx(1)=  -vb(2)*vc(3) + vb(3)*vc(2)
         dm_dx(2)=   vb(2)*vc(3) - vb(3)*vc(2) + va(2)*vc(3) - va(3)*vc(2)
         dm_dx(3)=   va(3)*vc(2) - va(2)*vc(3) + vb(2)*va(3) - vb(3)*va(2)
         dm_dx(4)=   va(2)*vb(3) - va(3)*vb(2)
       elseif(i==2)then
         dm_dx(1)=  -vb(3)*vc(1) + vb(1)*vc(3)
         dm_dx(2)=   vb(3)*vc(1) - vb(1)*vc(3) + va(3)*vc(1) - va(1)*vc(3)
         dm_dx(3)=  -va(3)*vc(1) + va(1)*vc(3) + vb(3)*va(1) - vb(1)*va(3)
         dm_dx(4)=  -vb(3)*va(1) + vb(1)*va(3)
       elseif(i==3)then
         dm_dx(1)=  -vb(1)*vc(2) + vb(2)*vc(1)
         dm_dx(2)=   vb(1)*vc(2) - vb(2)*vc(1) + va(1)*vc(2) - va(2)*vc(1)
         dm_dx(3)=  -va(1)*vc(2) + va(2)*vc(1) + va(2)*vb(1) - va(1)*vb(2)
         dm_dx(4)=  -va(2)*vb(1) + va(1)*vb(2)
       endif
       de_dx(1)=dm_dx(1)/sb; de_dx(2)=  m*(xyz(i,3)-xyz(i,2))/(sb2*sb)+dm_dx(2)/sb
       de_dx(4)=dm_dx(4)/sb; de_dx(3)= -m*(xyz(i,3)-xyz(i,2))/(sb2*sb)+dm_dx(3)/sb
       !!----------------------------------------------------------------------------
       dd_dx(1) = f*de_dx(1) - (e*(dh_dx(1) + j*l*dk_dx(1) + k*l*dj_dx(1) + j*k*dl_dx(1))/(g*g))
       dd_dx(2) = f*de_dx(2) - (e*(dh_dx(2) + j*l*dk_dx(2) + k*l*dj_dx(2) + j*k*dl_dx(2))/(g*g))
       dd_dx(3) = f*de_dx(3) - (e*(dh_dx(3) + j*l*dk_dx(3) + k*l*dj_dx(3) + j*k*dl_dx(3))/(g*g))
       dd_dx(4) = f*de_dx(4) - (e*(dh_dx(4) + j*l*dk_dx(4) + k*l*dj_dx(4) + j*k*dl_dx(4))/(g*g))
       grad(i,1) = dtheta_dd*dd_dx(1)
       grad(i,2) = dtheta_dd*dd_dx(2)
       grad(i,3) = dtheta_dd*dd_dx(3)
       grad(i,4) = dtheta_dd*dd_dx(4)
       !rewriting this part below
       do jitr = 1, 4
        if(isnan(grad(i,jitr))) nan = .true.
       enddo
       if(nan) then
       !if(any(isnan(grad(i,1:4))))then
           write(55,*)dtheta_dd,theta,d,e,g,dd_dx(1:4)
           stop
       endif
   enddo
   !--------------------------------------
  end subroutine dihedral_gradient
!----------------------------------------------------------------------------------------------------------------------------
  subroutine initialize_velocity(model,vel,nres,t)
   !------------------------------------
   implicit none
   !-------------------------------------
   integer,intent(in)::nres
   real,intent(in)::t
   real,dimension(:,:),pointer::model
   real,dimension(:,:),pointer::vel
   !-------------------------------------
   integer::k
   integer::i
   real::lamba
   real::curr_t
   real::ke
   real,dimension(3)::vcom
   !-------------------------------------
   call random_number(vel)
   vel = vel - 0.5
   call remove_net_momentum(model,velocity,nres)
   !--------------------------------------
  end subroutine
!----------------------------------------------------------------------------------------------------------------------------
  subroutine remove_net_momentum(xyz,vel,nres)
   !--------------------------------------
   implicit none
   !--------------------------------------
   integer,intent(in)::nres
   real,dimension(:,:),pointer::vel
   real,dimension(:,:),intent(in)::xyz
   !--------------------------------------
   integer::ierr
   integer::ires
   real::ixx
   real::iyy
   real::izz
   real::ixy
   real::ixz
   real::iyz
   real,dimension(3)::r
   real,dimension(3)::com
   real,dimension(3)::angular
   real,dimension(3)::linear
   real,dimension(3)::omega
   real,dimension(3)::wxr
   real,dimension(3)::diag
   real,dimension(3)::l
   real,dimension(3,3)::itensor
   !--------------------------------------
   com=0
   angular=0
   linear=0
   do ires=1,nres
     com = com + xyz(:,ires)
     linear = linear + vel(:,ires)
   enddo
   com = com/real(nres)
   linear = linear/real(nres)
   do ires=1,nres
     r = xyz(:,ires) - com
     call crossprod(r,vel(:,ires),l)
     angular = angular + l
   enddo 
   ixx=0; iyy=0; izz=0; ixy=0; ixz=0; iyz=0
   do ires=1,nres
     r = xyz(:,ires) - com
     ixx = ixx + r(2)*r(2) + r(3)*r(3)
     iyy = iyy + r(1)*r(1) + r(3)*r(3)
     izz = izz + r(1)*r(1) + r(2)*r(2)
     ixy = ixy - r(1)*r(2)
     ixz = ixz - r(1)*r(3)
     iyz = iyz - r(2)*r(3)   
   enddo 
   itensor(1,1)=ixx;  itensor(2,2)=iyy;  itensor(3,3)=izz
   itensor(1,2)=ixy;  itensor(2,1)=ixy
   itensor(1,3)=ixz;  itensor(3,1)=ixz
   itensor(2,3)=iyz;  itensor(3,2)=iyz
   call choldc(itensor,diag,3,ierr)
   call cholsl(itensor,3,diag,angular,omega)
   if(ierr/=0) return
   do ires=1,nres
     r = xyz(:,ires) - com
     call crossprod(omega,r,wxr)
     vel(:,ires) = vel(:,ires)-wxr-linear 
   enddo
   !------------------------------------------
  end subroutine remove_net_momentum
!----------------------------------------------------------------------------------------------------------------------------
  subroutine choldc(a,p,n,izero)
   !------------------------------------------
   implicit none
   !------------------------------------------
   !* A = matrix  N = size  P = diagonal
   !------------------------------------------
   integer,intent(in)::n
   integer,intent(out)::izero
   real,dimension(n),intent(out)::p
   real,dimension(n,n),intent(inout)::a
   !------------------------------------------
   integer::i
   integer::j
   integer::k
   real::ssum
   !------------------------------------------
   izero = 0
   do i=1,n
     do j=i,n
       ssum = a(j,i)
       do k=i-1,1,-1
         ssum = ssum - a(k,i)*a(k,j)
       enddo
       if (i.eq.j) then
         if (ssum.le.0.0000) then
           if (ssum.lt.0.0000) then
             write(*,*) ' subroutine choldc :: failed on row= ',i
             write(*,*) '  Sum= ',ssum
             write(*,*) '  A(i,i) before summing= ',a(j,i)
           endif
           izero = i
           return 
         endif
         p(i) = sqrt(ssum)
       else
         a(i,j) = ssum/p(i)
       endif
     enddo
   enddo
   !------------------------------------------
  end subroutine choldc
!----------------------------------------------------------------------------------------------------------------------------
  subroutine cholsl(a,n,p,b,x)
   !-----------------------------------------------
   implicit none
   !-----------------------------------------------
   !* A is from subroutine choldc  N is size  P is from CHOLDC
   !* B is known, X is the output vector 
   !-----------------------------------------------
   integer,intent(in)::n
   real,dimension(n),intent(in)::p
   real,dimension(n),intent(in)::b
   real,dimension(n),intent(out)::x
   real,dimension(n,n),intent(in)::a
   !-----------------------------------------------
   integer::i
   integer::j
   real::ssum
   !-----------------------------------------------
   do i=1,n
     ssum = b(i)
     do j=i-1,1,-1
       ssum = ssum - a(j,i)*x(j)
     enddo
     x(i) = ssum/p(i)
   enddo
   do i=n,1,-1
     ssum = x(i)
     do j=i+1,n
       ssum = ssum - a(i,j)*x(j)
     enddo
     x(i) = ssum/p(i)
   enddo
   !------------------------------------------
  end subroutine cholsl
!----------------------------------------------------------------------------------------------------------------------------
  subroutine crossprod(v1,v2,v3)
   !------------------------------------------------
   implicit none
   !------------------------------------------------
   real,intent(in)::v1(3)
   real,intent(in)::v2(3)
   real,intent(out)::v3(3)
   !------------------------------------------------
   v3(1) = v1(2)*v2(3) - v1(3)*v2(2)
   v3(2) = v1(3)*v2(1) - v1(1)*v2(3)
   v3(3) = v1(1)*v2(2) - v1(2)*v2(1)
   !------------------------------------------
  end subroutine crossprod
!----------------------------------------------------------------------------------------------------------------------------
  subroutine velocity_verlet_1(acc,xyz,vel,nres,dt)
   !------------------------------------------
   implicit none
   !------------------------------------------
   integer,intent(in)::nres
   real,intent(in)::dt
   real,dimension(:,:),pointer::acc
   real,dimension(:,:),intent(inout)::xyz
   real,dimension(:,:),pointer::vel
   !------------------------------------------
   integer::ires
   real::dt2
   real::dtsq2
   !------------------------------------------
   !half_dt=dt*0.5
   !do ires=1, nres
   !  vel(:,ires) = vel(:,ires) + half_dt*acc(:,ires)   
   !  xyz(:,ires) = xyz(:,ires) + dt*vel(:,ires)
   !enddo
   !------------------------------------------
   dt2=dt*0.5
   dtsq2=dt2*dt
   do ires=1, nres
     !if(ires/=go_s .and. ires/=go_e) then
     xyz(:,ires) = xyz(:,ires) + dt*vel(:,ires) + dtsq2*acc(:,ires)
     vel(:,ires) = vel(:,ires)                  +   dt2*acc(:,ires) !I am almost postive this should be dt not dt2
     !endif  
   enddo
   !------------------------------------------
  end subroutine velocity_verlet_1
!----------------------------------------------------------------------------------------------------------------------------
  subroutine velocity_verlet_2(acc,vel,nres,dt) !adjusts the current velocity for the current acceleration 
   implicit none
   !------------------------------------------
   integer,intent(in)::nres
   real,intent(in)::dt
   real,dimension(:,:),pointer::vel
   real,dimension(:,:),pointer::acc
   !------------------------------------------
   integer::ires
   real::dt2
   !------------------------------------------
   !half_dt=dt*0.5
   !do ires=1,nres
   !  vel(:,ires) = vel(:,ires) + half_dt*acc(:,ires)   
   !enddo
   !------------------------------------------
   dt2=dt*0.5
   do ires=1,nres
     if(ires/=go_s .and. ires/=go_e) then
     vel(:,ires) = vel(:,ires) + dt2*acc(:,ires) !same here dt2 should really be dt
     endif  
   enddo
   !------------------------------------------
  end subroutine velocity_verlet_2
!----------------------------------------------------------------------------------------------------------------------------
  subroutine set_target_dihedral(dih,nres)
   !------------------------------------------
   implicit none
   !------------------------------------------
   integer,intent(in)::nres
   real,dimension(:),pointer::dih
   !------------------------------------------
   integer::ires
   !------------------------------------------
   do ires=1,nres
     if(dih(ires)< 180.0) then
       dih(ires) = dih(ires) + 180.0
     else
       dih(ires) = dih(ires) - 180.0
     endif
   enddo
   !------------------------------------------
  end subroutine set_target_dihedral
!----------------------------------------------------------------------------------------------------------------------------
!  subroutine check_torq(xyz,fshift,natm)
!   !--------------------------------------
!   implicit none
!   !--------------------------------------
!   integer,intent(in)::natm
!   real,dimension(3,natm),intent(in)::xyz
!   real,dimension(3,natm),intent(in)::fshift
!   !--------------------------------------
!   integer::iatm
!   real,dimension(3)::com
!   real,dimension(3)::linear
!   real,dimension(3)::torq
!   real,dimension(3)::rvec
!   real,dimension(3)::tvec
!   !--------------------------------------
!   com=0
!   torq=0
!   linear=0
!   do iatm=1,natm
!     com    =  com    + xyz(:,iatm)
!     linear =  linear + fshift(:,iatm)
!   enddo
!   com = com /real(natm)
!   do iatm=1,natm
!     rvec = xyz(:,iatm)  - com
!     call crosprod(rvec,fshift(:,iatm),tvec)
!     torq = torq + tvec
!   enddo 
!   write(sample_unit,'(a,6f10.5)')"=====> ",linear,torq
!   !--------------------------------------
!  end subroutine check_torq
!----------------------------------------------------------------------------------------------------------------------------
end program cago
