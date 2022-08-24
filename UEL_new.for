       module kvisual
       implicit none
       real*8 sigout(10000,4,4)
       save
       end module    

	   SUBROUTINE UEL(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,
     1 PROPS,NPROPS,COORDS,MCRD,NNODE,U,DU,V,A,JTYPE,TIME,DTIME,
     2 KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,PREDEF,NPREDF,
     3 LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,NJPROP,PERIOD)
C
		use kvisual
      INCLUDE 'ABA_PARAM.INC'
C
	  
       dimension RHS(MLVARX,*),AMATRX(NDOFEL,NDOFEL),PROPS(*),
     1 SVARS(NSVARS),ENERGY(8),COORDS(MCRD,NNODE),U(NDOFEL),
     2 DU(MLVARX,*),V(NDOFEL),A(NDOFEL),TIME(2),PARAMS(*),
     3 JDLTYP(MDLOAD,*),ADLMAG(MDLOAD,*),DDLMAG(MDLOAD,*),
     4 PREDEF(2,NPREDF,NNODE),LFLAGS(*),JPROPS(*)

	   PARAMETER (maxDOF=2 , ndof=2 , ndim=2 ,ninpt=4,nelem=16,nsdv=4,nsvint=4,ntens=4)
       dimension dshape(ndim,4),xjac(maxDOF,maxDOF),xjaci(maxDOF,maxDOF),
     1 bmat(3,2*nnode),bmatT(2*nnode,3),kmat(2*nnode,2*nnode),C(3,3),
     2 deps(3),stress(3),SIG(4),SDV(4),statevLocal(nsvint)
	   dimension db(3,2*nnode),XYDerivation(ndim,4),force(8)

		rhs(:, 1)=0.d0
		amatrx=0.d0
		C=0.d0		
				
		
		force=0.d0
		
		!Elasticity matrix
		E = props(1)
		xnu = props(2)
		xa = E/(1.d0- xnu**2)
		xb = xnu
		xc = (1.d0-xnu)/2.d0
		
		  C(1, 1) = xa
		  C(1, 2) = xa*xb
		  C(2, 1) = C(1, 2)
		  C(2, 2) = C(1, 1)
		  C(3, 3) = xa*xc

        do kintk = 1, ninpt	 
		deps=0.	
		statevLocal=0.d0
        call shapeFunc(kintk,ninpt,nnode,ndim,shapef,dshape)
		
	    djac = 1.d0
        call jacobian(ndim,COORDS,nnode,djac,dshape,xjaci)
	    
        call Bmatrix(xjaci,dshape,nnode,ndim,bmat,bmatT)

	    call statevar(kintk,nsvint,SVARS,NSVARS,statevLocal,1)

		SIG=0.d0
		stress=0.d0
		do k1=1,nsdv
			SIG(k1) = statevLocal(k1)
        end do
		
	    stress(1)=SIG(1)
		stress(2)=SIG(2)
		stress(3)=SIG(4)
		   
		call straininc(ntens,ndof,ndim,nnode,MLVARX,bmat,DU,deps)		


		stress=stress+matmul(C,deps)
		force=0.d0
		force=matmul(bmatT,stress)
		force=force*djac
		do k1=1,8
			
				RHS(k1, 1) = RHS(k1, 1) - force(k1)
			
		end do
		
		db=0.
	   	do i=1,3
			do j=1,8
				do k=1,3
					db(i,j)=db(i,j)+C(i,k)*bmat(k,j)*djac
				end do	
			end do
		end do
		
	   
		do i=1,8
			do j=1,8
				do k=1,3
					AMATRX(i,j)=AMATRX(i,j)+bmatT(i,k)*db(k,j)
				end do	
			end do
		end do
		
		SIG(1)=stress(1)
		SIG(2)=stress(2)
		SIG(3)=0.d0
		SIG(4)=stress(3)
		
	    do k1=1,nsdv
			statevLocal(k1) = SIG(k1)
		end do

        call statevar(kintk,nsvint,SVARS,NSVARS,statevLocal,0)
	   
	

		do k1=1,nsdv
			sigout(jelem,kintk,k1)=SIG(k1)
		end do
		
		end do       


      RETURN
      END
	  
	  
	  
******************************************************
** shape functions                                   *
******************************************************	
		subroutine shapeFunc(kintk,ninpt,nnode,ndim,shapef,dshape)
		
		include 'aba_param.inc'
		
		parameter (gaussCoord=0.577350269189626)
		dimension shapef(4),dshape(ndim,4),coord24(2,4)
		
		data coord24  /1. , 1. ,
	2				  -1. , 1. ,
	3				   1. , -1.,
	4				   -1. , -1./
			
		shapef=0.
		dshape=0.
		r=coord24(1,kintk)*gaussCoord
		s=coord24(2,kintk)*gaussCoord
		
		! shape functions
		
		shapef(1)=0.25*(1.+r)*(1.+s)
		shapef(2)=0.25*(1.-r)*(1.+s)
        shapef(3)=0.25*(1.-r)*(1.-s)
        shapef(4)=0.25*(1+r)*(1.-s)
	
		! derivation to r
		
		dshape(1,1)=0.25*(1.+s)
		dshape(1,2)=-0.25*(1.+s)
        dshape(1,3)=-0.25*(1.-s)
        dshape(1,4)=0.25*(1.-s)	
	
	
	    ! derivation to s
		
		dshape(2,1)=0.25*(1.+r)
		dshape(2,2)=0.25*(1.-r)
        dshape(2,3)=-0.25*(1.-r)
        dshape(2,4)=-0.25*(1.+r)
	
	
		return
		end
		
******************************************************
** Jacobian                                          *
******************************************************
		subroutine jacobian(ndim,coords,nnode,djac,dshape,xjaci)
		
		include 'aba_param.inc'
		
		dimension xjac(ndim,ndim),xjaci(ndim,ndim)
		dimension coords(2,nnode),dshape(ndim,nnode)

		xjac=0.d0
		xjaci=0.d0

		
        do inod= 1, nnode
          do kdim = 1, ndim
            do jdim = 1, ndim
              xjac(jdim,kdim) = xjac(jdim,kdim) + 
     1        dshape(jdim,inod)*coords(kdim,inod)      
            end do
          end do 
        end do
		
		
		djac=xjac(1,1)*xjac(2,2)-xjac(1,2)*xjac(2,1)
		 
		if (djac .gt. 0.) then
			xjaci(1,1)=xjac(2,2)/djac
			xjaci(2,2)=xjac(1,1)/djac
			xjaci(1,2)=-xjac(1,2)/djac
			xjaci(2,1)=-xjac(2,1)/djac			
		end if
		
		return
		end
******************************************************
** B-matrix                                          *
******************************************************
		subroutine Bmatrix(xjaci,dshape,nnode,ndim,bmat,bmatT)
		
		include 'aba_param.inc'
		
		dimension xjaci(ndim,ndim),dshape(ndim,4),XYDerivation(ndim,4)
		dimension bmat(3,2*nnode),bmatT(2*nnode,3)

				XYDerivation=0.d0

				bmat=0.d0
				bmatT=0.d0


					XYDerivation=matmul(xjaci,dshape)

				

		bmat(1,1)=XYDerivation(1,1)
		bmat(1,3)=XYDerivation(1,2)
		bmat(1,5)=XYDerivation(1,3)
		bmat(1,7)=XYDerivation(1,4)
		
		bmat(2,2)=XYDerivation(2,1)
		bmat(2,4)=XYDerivation(2,2)
		bmat(2,6)=XYDerivation(2,3)
		bmat(2,8)=XYDerivation(2,4)
		
		bmat(3,1)=XYDerivation(2,1)
		bmat(3,2)=XYDerivation(1,1)
		bmat(3,3)=XYDerivation(2,2)
		bmat(3,4)=XYDerivation(1,2)
		bmat(3,5)=XYDerivation(2,3)
		bmat(3,6)=XYDerivation(1,3)
		bmat(3,7)=XYDerivation(2,4)
		bmat(3,8)=XYDerivation(1,4)
		
		do i=1,3
			do j=1,2*nnode
				bmatT(j,i)=bmat(i,j)
			end do
		end do
		
		return
		end
		
c*****************************************************************
		  subroutine straininc(ntens,ndof,ndim,nnode,mlvarx,bmat,du,deps)

			include 'aba_param.inc'
		  dimension deps(3),bmat(3,2*nnode),du(mlvarx,1),xdu(2) 
		  
			
		  deps = 0.d0
		  xdu=0.d0
			
		  do nodi = 1, nnode
			   
		   incr_row = (nodi - 1)*ndof
			
		   do i = 1, ndof         
				xdu(i)= du(i + incr_row,1)
		   end do
			dNidx = bmat(1,1 + incr_row)
			dNidy = bmat(2,1 + incr_row + 1)
			
		   deps(1) = deps(1) + dNidx*xdu(1)
		   deps(2) = deps(2) + dNidy*xdu(2)
		   deps(3) = deps(3) + dNidy*xdu(1) + dNidx*xdu(2)
					   
		  end do

		  return
		  end


		
************************************************************************		
		 subroutine statevar(npt,nsvint,statev,NSVARS,statev_ip,icopy)
	

		  include 'aba_param.inc'

		  dimension statev(NSVARS),statev_ip(nsvint)

		  isvinc= (npt-1)*nsvint     

		  if (icopy .eq. 1) then

			do i = 1, nsvint
			  statev_ip(i)=statev(i+isvinc)
			end do
c
		  else
c
			do i = 1, nsvint
			  statev(i+isvinc)=statev_ip(i)
			end do
		  end if

		  return
		  end

************************************************************************
      subroutine umat(stress,statev,ddsdde,sse,spd,scd,
     1 rpl,ddsddt,drplde,drpldt,
     2 stran,dstran,time,dtime,temp,dtemp,predef,dpred,cmname,
     3 ndi,nshr,ntens,nstatv,props,nprops,coords,drot,pnewdt,
     4 celent,dfgrd0,dfgrd1,noel,npt,layer,kspt,kstep,kinc)
C
		use kvisual
      include 'aba_param.inc'
C
      character*80 cmname
      dimension stress(ntens),statev(nstatv),
     1 ddsdde(ntens,ntens),
     2 ddsddt(ntens),drplde(ntens),
     3 stran(ntens),dstran(ntens),time(2),predef(1),dpred(1),
     4 props(nprops),coords(3),drot(3,3),dfgrd0(3,3),dfgrd1(3,3)
	  PARAMETER (maxDOF=2 , ndof=2 , ndim=2 ,ninpt=4,nelem=16,nsdv=4)
	  integer kelem
      


		ddsdde=0.d0

      kelem=int(noel-nelem)
	 
      do k1 = 1, nsdv
        statev(k1) = sigout(kelem,npt,k1)
      end do 

      return
      end 