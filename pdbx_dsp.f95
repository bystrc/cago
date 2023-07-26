module pdbx
  !-------------------------------
  implicit none
  !-------------------------------
  type,public:: residue
    real,dimension(3,22)::xyz
  end type
  !-------------------------------
  private
  integer,dimension(21)::atom_number=(/5,6,8,9,11,4,10,8,9,8,8,8,7,9,11,6,7,7,14,12,22/)
  character(len=3),dimension(22)::three=(/'ALA','CYS','ASP','GLU','PHE','GLY','HIS','ILE',&
     'LYS','LEU','MET','ASN','PRO','GLN','ARG','SER','THR','VAL',&
     'TRP','TYR','CRO','XXX'/)
  character(len=1),dimension(22)::res1=(/'A','C','D','E','F','G','H','I','K','L','M','N',&
     'P','Q','R','S','T','V','W','Y','Z','X'/)
  character(4), dimension(28,21)::atype=reshape(&
   (/' N  ',' CA ',' C  ',' O  ',' CB ','    ','    ','    ','    ','    ','    ','    ','    ','    ',& ! A
     '    ','    ','    ','    ','    ','    ','    ','    ','    ','    ','    ','    ','    ','    ',& ! A
     ' N  ',' CA ',' C  ',' O  ',' CB ',' SG ','    ','    ','    ','    ','    ','    ','    ','    ',& ! C
     '    ','    ','    ','    ','    ','    ','    ','    ','    ','    ','    ','    ','    ','    ',& ! C
     ' N  ',' CA ',' C  ',' O  ',' CB ',' CG ',' OD1',' OD2','    ','    ','    ','    ','    ','    ',& ! D
     '    ','    ','    ','    ','    ','    ','    ','    ','    ','    ','    ','    ','    ','    ',& ! D
     ' N  ',' CA ',' C  ',' O  ',' CB ',' CG ',' CD ',' OE1',' OE2','    ','    ','    ','    ','    ',& ! E
     '    ','    ','    ','    ','    ','    ','    ','    ','    ','    ','    ','    ','    ','    ',& ! E
     ' N  ',' CA ',' C  ',' O  ',' CB ',' CG ',' CD1',' CD2',' CE1',' CE2',' CZ ','    ','    ','    ',& ! F
     '    ','    ','    ','    ','    ','    ','    ','    ','    ','    ','    ','    ','    ','    ',& ! F
     ' N  ',' CA ',' C  ',' O  ','    ','    ','    ','    ','    ','    ','    ','    ','    ','    ',& ! G
     '    ','    ','    ','    ','    ','    ','    ','    ','    ','    ','    ','    ','    ','    ',& ! G
     ' N  ',' CA ',' C  ',' O  ',' CB ',' CG ',' ND1',' CD2',' CE1',' NE2','    ','    ','    ','    ',& ! H
     '    ','    ','    ','    ','    ','    ','    ','    ','    ','    ','    ','    ','    ','    ',& ! H
     ' N  ',' CA ',' C  ',' O  ',' CB ',' CG1',' CG2',' CD1','    ','    ','    ','    ','    ','    ',& ! I
     '    ','    ','    ','    ','    ','    ','    ','    ','    ','    ','    ','    ','    ','    ',& ! I
     ' N  ',' CA ',' C  ',' O  ',' CB ',' CG ',' CD ',' CE ',' NZ ','    ','    ','    ','    ','    ',& ! K
     '    ','    ','    ','    ','    ','    ','    ','    ','    ','    ','    ','    ','    ','    ',& ! K
     ' N  ',' CA ',' C  ',' O  ',' CB ',' CG ',' CD1',' CD2','    ','    ','    ','    ','    ','    ',& ! L
     '    ','    ','    ','    ','    ','    ','    ','    ','    ','    ','    ','    ','    ','    ',& ! L
     ' N  ',' CA ',' C  ',' O  ',' CB ',' CG ',' SD ',' CE ','    ','    ','    ','    ','    ','    ',& ! M
     '    ','    ','    ','    ','    ','    ','    ','    ','    ','    ','    ','    ','    ','    ',& ! M
     ' N  ',' CA ',' C  ',' O  ',' CB ',' CG ',' OD1',' ND2','    ','    ','    ','    ','    ','    ',& ! N 
     '    ','    ','    ','    ','    ','    ','    ','    ','    ','    ','    ','    ','    ','    ',& ! N
     ' N  ',' CA ',' C  ',' O  ',' CB ',' CG ',' CD ','    ','    ','    ','    ','    ','    ','    ',& ! P
     '    ','    ','    ','    ','    ','    ','    ','    ','    ','    ','    ','    ','    ','    ',& ! P
     ' N  ',' CA ',' C  ',' O  ',' CB ',' CG ',' CD ',' OE1',' NE2','    ','    ','    ','    ','    ',& ! Q             
     '    ','    ','    ','    ','    ','    ','    ','    ','    ','    ','    ','    ','    ','    ',& ! Q
     ' N  ',' CA ',' C  ',' O  ',' CB ',' CG ',' CD ',' NE ',' CZ ',' NH1',' NH2','    ','    ','    ',& ! R
     '    ','    ','    ','    ','    ','    ','    ','    ','    ','    ','    ','    ','    ','    ',& ! R
     ' N  ',' CA ',' C  ',' O  ',' CB ',' OG ','    ','    ','    ','    ','    ','    ','    ','    ',& ! S
     '    ','    ','    ','    ','    ','    ','    ','    ','    ','    ','    ','    ','    ','    ',& ! S
     ' N  ',' CA ',' C  ',' O  ',' CB ',' OG1',' CG2','    ','    ','    ','    ','    ','    ','    ',& ! T
     '    ','    ','    ','    ','    ','    ','    ','    ','    ','    ','    ','    ','    ','    ',& ! T
     ' N  ',' CA ',' C  ',' O  ',' CB ',' CG1',' CG2','    ','    ','    ','    ','    ','    ','    ',& ! V
     '    ','    ','    ','    ','    ','    ','    ','    ','    ','    ','    ','    ','    ','    ',& ! V
     ' N  ',' CA ',' C  ',' O  ',' CB ',' CG ',' CD1',' CD2',' NE1',' CE2',' CE3',' CZ2',' CZ3',' CH2',& ! W
     '    ','    ','    ','    ','    ','    ','    ','    ','    ','    ','    ','    ','    ','    ',& ! V
     ' N  ',' CA ',' C  ',' O  ',' CB ',' CG ',' CD1',' CD2',' CE1',' CE2',' CZ ',' OH ','    ','    ',& ! Y
     '    ','    ','    ','    ','    ','    ','    ','    ','    ','    ','    ','    ','    ','    ',& ! V
     ' N1 ',' CA1',' CB1',' CG1',' OG1',' C1 ',' N2 ',' N3 ',' C2 ',' O2 ',' CA2',' CA3',' C3 ',' O3 ',& ! Z
     ' CB2',' CG2',' CD1',' CD2',' CE1',' CE2',' CZ ',' OH ','    ','    ','    ','    ','    ','    ' & ! Z
     /),(/28,21/))
  public::pdbx_read_protein,pdbx_write_protein
  !-------------------------------------------------------------------------------
 contains
!----------------------------------------------------------------------------------------------------------------------------
  subroutine pdbx_read_protein(pdbfile,chain,protein,dsbonds,ndsbonds,pdbseq,calpha,nres)!,start,end)
   !-------------------------------
   implicit none
   !-------------------------------
   type(residue),dimension(:),pointer::protein
   integer,dimension(:,:),pointer::dsbonds
   integer,intent(out)::ndsbonds
   character(len=*),intent(in)::pdbfile
   character(len=1),intent(inout)::chain
   integer,intent(inout)::nres
   !integer,intent(in):: start
   !integer,intent(in)::end
   integer,dimension(:),intent(inout)::pdbseq
   real,dimension(:,:),pointer::calpha

   !-------------------------------
   character(len=1)::altloc
   character(len=100)::aline
   character(len=6)::last
   character(len=10)::curres
   character(len=10)::qres
   integer::i
   integer::j
   integer::k
   integer::n
   !integer::s_res
  ! integer::e_res
   integer::dsresidue
   integer::iatom
   integer::jarg
   integer::ios
   integer::ires
   integer::idsbond
   integer::dsbond_n
   integer::curratom
   integer::iline
   integer::atype
   integer::restype
   !-------------------------------------------------------------------------------
   open(2,file=pdbfile,iostat=ios)
   if(ios/=0)then
     write(*,*) "Cannot open pdbfile ",trim(pdbfile),"!"
     stop
   endif
   !-------------------------------------------------------------------------------
   ires=0;curres=" "; idsbond = 0;
   !-------------------------------------------------------------------------------
   do 
     read(2,'(a)',iostat=ios) aline
     !write(*,*)aline,'  ',ios
     if(ios/=0) exit
     if(aline(1:6)=='ENDMDL')exit              !! use first model of NMR structures!
     if(aline(1:3)=="TER".and.ires/=0)exit     !! Atoms with the correct chain id but not part of the chain! 
     if(aline(1:6) == 'SSBOND')then
     !read(aline(18:21),'(I4)') s_res
     !read(aline(32:35),'(I4)') e_res 
      !if((s_res>=start .and. s_res<=end) &
!&     .and. (e_res>=start .and. e_res<=end)) then
     idsbond = idsbond + 1
!      endif
     endif
     if(aline(1:6)/='ATOM  '.and.aline(1:6)/='HETATM') cycle
     if(aline(22:22)/=chain.and.aline(22:22)/=" ")cycle
     if(aline(22:22)==" ")chain="_"
     if(aline(17:17)/= "A" .and. aline(17:17)/=" ") cycle !only take the first position
     if(aline(13:16)/=' CA '.and.aline(18:20)/='CRO')cycle
     if(aline(18:20)=='CRO'.and.aline(13:16)/=' CA1')cycle
     if(aline(18:27)==curres)cycle             !! Only one CA per residue allowed!
    ! if(ires>end .and. end>0) exit
     curres=aline(18:27)
     ires = ires + 1
     !if(ires<start) cycle ! read only after start
     !write(*,*)ires," ",curres
     restype = 22
     do j=1,21
       if(aline(18:20)==three(j)) then
         restype=j; exit
       endif
     enddo
     pdbseq(ires) = j
     if(j==22) write(*,*) " subroutine pdbx_read_protein :: unrecognized residue type ",trim(aline(18:20))
   
     !write(*,*) ires,'*'
  enddo
   !-------------------------------------------------------------------------------
!   nres=ires-start
    nres=ires
   ndsbonds= idsbond
   !write(*,*)trim(pdbfile)," ====> ",nres," residues."
   !-------------------------------------------------------------------------------
   if(nres==0) then
     write(*,*) "Parsing error, no residues read!"
     stop
   endif
   !write(*,*) ndsbonds, "disulfide bonds found"
   !write(*,*) nres, "residues selected"
   !-------------------------------------------------------------------------------
   allocate(protein(nres))
   allocate(dsbonds(2,ndsbonds))
   !-------------------------------------------------------------------------------
   do ires=1,nres
     protein(ires)%xyz=999.
   enddo
   rewind(2)
   !-------------------------------------------------------------------------------
   curres=" "
   ires=0 
   idsbond=1
   !-------------------------------------------------------------------------------
   do 
     iatom=0
     read(2,'(a)',iostat=ios) aline
     if(ios/=0) exit
     if(aline(1:6)=='ENDMDL')exit  
     if(aline(1:3)=="TER".and.ires/=0)exit                   
     if(aline(1:6)=="SSBOND" .and. (aline(12:14) /= "CYS" .or. aline(26:28)/= "CYS")) then 
       write(*,*) "only cysteine forms disulfide bonds! "
       cycle
     endif
     if(aline(1:6)=="SSBOND" .and. aline(12:14)=="CYS") then
        read (aline(18:21),'(I4)') dsresidue
        dsbonds(1,idsbond) = dsresidue
        read (aline(32:35), '(I4)') dsresidue
        dsbonds(2,idsbond) = dsresidue
        idsbond=idsbond+1
     endif 
     if(aline(1:6)/='ATOM  '.and.aline(1:6)/='HETATM') cycle
     if(aline(22:22)/=chain.and.aline(22:22)/=" ")cycle
     if(aline(18:27)==curres)cycle
    ! if(ires<start) cycle ! read only after start
    ! if(ires>end) exit

     !-----------------------------------------------------------------------------
     restype = 22
     do j=1,21
       if(aline(18:20)==three(j)) then
          restype=j; exit
       endif
     enddo
     if(restype==22) then
        write(*,*) "Cannot parse unrecognized residue type!",curres," ",aline(18:26)," ",trim(pdbfile)," ",chain
        stop
     endif
     !-----------------------------------------------------------------------------
     atype = idatom(aline(13:16),aline(18:20))
     !-----------------------------------------------------------------------------
     if(atype==-1)cycle                                      
     !-----------------------------------------------------------------------------
     ires = ires + 1; iatom = iatom + 1        
     !-----------------------------------------------------------------------------
     !if(ires>e_res)then
     if(ires>nres) then
        write(*,*)"Missing atoms?(1) ",curres," ",aline(18:26)," ",trim(pdbfile)," ",chain
        stop
     endif
     !-----------------------------------------------------------------------------
     read(aline(31:54),'(3f8.3))')protein(ires)%xyz(1:3,atype) 
     !-----------------------------------------------------------------------------
     curres=aline(18:27)                                     
     !-----------------------------------------------------------------------------
     do while(iatom/=atom_number(restype))   
       !-----------------------------------------------------------------------------
       altloc=" "
       !-----------------------------------------------------------------------------
       read(2,'(a)',iostat=ios)aline
       !-----------------------------------------------------------------------------
       if(aline(1:6)/='ATOM  '.and.aline(1:6)/='HETATM') cycle
       !-----------------------------------------------------------------------------
       qres=aline(18:27)                                     
       !-----------------------------------------------------------------------------
       if(curres/=qres)then
         write(*,*)"Missing atoms?(2) ",curres," ",trim(pdbfile)," ",chain
         backspace(2)
         exit
       endif
       !-----------------------------------------------------------------------------
       atype = idatom(aline(13:16),aline(18:20))
       !-----------------------------------------------------------------------------
       if(atype==-1)cycle
       !-----------------------------------------------------------------------------
       altloc=aline(17:17)                              
       !-----------------------------------------------------------------------------
       if(all(protein(ires)%xyz(1:3,atype)==999.))then
         read(aline(31:54),'(3f8.3))')protein(ires)%xyz(1:3,atype) 
         iatom = iatom + 1
       elseif(altloc/=" ")then
         cycle                                           
       else
         write(*,*)"Attempting to overwrite atomic coordinates!",ires,atype,trim(pdbfile)," ,",chain
         write(*,*)aline(18:26)
         stop
       endif
       !-----------------------------------------------------------------------------
     enddo
   enddo
   close(2)
   !-------------------------------------------------------------------------------
   allocate(calpha(3,nres))
   !-------------------------------------------------------------------------------
   do ires=1,nres
     calpha(:,ires)=protein(ires)%xyz(:,2)
   enddo
   !-------------------------------------------------------------------------------
  end subroutine
!----------------------------------------------------------------------------------------------------------------------------
  integer function idatom(atom,res)
   !-------------------------------
   implicit none
   !-------------------------------
   character(len=4),intent(in)::atom
   character(len=3),intent(in)::res
   !-------------------------------
   integer::i
   integer::j
   integer::k
   !-------------------------------
   do j=1,21
     if(res==three(j)) then 
       k=j;exit
     endif
   enddo
   do i=1,28
     if(atom==atype(i,k))exit
   enddo
   if(i==29)then
     idatom = -1
   else
     idatom=i
   endif
   !-------------------------------
  end function
!----------------------------------------------------------------------------------------------------------------------------
  subroutine pdbx_write_protein(iunit,protein,nres,seq,chain)
   !-------------------------------
   implicit none
   !-------------------------------
   type(residue),dimension(:),intent(in)::protein
   character(len=1),intent(in)::chain
   integer,intent(in)::nres
   integer,intent(in)::iunit
   integer,dimension(:),intent(in)::seq
   !-------------------------------
   integer::ires
   integer::iatom
   integer::j
   integer::jatom
   integer::restype
   !-------------------------------
   jatom=0
   do ires=1,nres
     do iatom=1,atom_number(seq(ires))
       if(protein(ires)%xyz(1,iatom)/=999.) then
         jatom=jatom+1
         write(iunit,'("ATOM",3X,i4,1X,a4,1x,a3,1x,a1,i4,4X,3F8.3)')&
         jatom,atype(iatom,seq(ires)),three(seq(ires)),chain,ires,protein(ires)%xyz(:,iatom)
       endif
     enddo
   enddo
  !-------------------------------
  end subroutine pdbx_write_protein
!----------------------------------------------------------------------------------------------------------------------------
end module pdbx
