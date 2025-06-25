subroutine thetaface(ux,uy,uz,x1,y1,z1,tan_th2_1,tan_th2_2)

  use closest_wall
  !
  ! Computes the shortest distance (t) to constant theta (polar angle)
  ! surface set of (tan(theta1),tan(theta2)) from the point (x1,y1,z1)
  ! along the line defined by the directional cosines (ux,uy,uz).
  ! Intersections with negative "t" are NOT selected (i.e., are
  ! OPPOSITE of desired direction of propagation).
  !
  ! NOTE:  constant theta surfaces are assumed to be centered at origin
  ! with z-axis defining theta=0.  surfaces are cylindrical cones.
  ! (x/a)**2 + (y/a)**2 = (z/c)**2, tan(theta) = a/c
  ! ->  x**2 + y**2 - (z*tan(theta))**2 = 0.
  !
  !
  ! Output values:
  ! iface = which constant radius surface is closest
  ! = 0 for theta1, 1 for theta2, and -1 for NEITHER
  ! ifound = flag for non-negative
  ! = .T. for positive "t" found, .F. NOT.
  ! t = smallest distance for surface intersection
  ! (negative value indicated desired interesction does NOT exist...
  ! but user should test for this case using "ifound".
  !
  ! Input values:
  ! ux,uy,uz = directional cosines along x,y,z directions
  ! x1,y1,z1 = point from which distance is calculated
  ! tan_th2_1,_2 = (tan(theta))**2 for two constant surfaces
  ! NOTE: ***** it is necessary to flag tan_th2 array using a -1.
  ! to indicate theta=90.
  ! Program ASSUMES that theta=90 can occur in ONLY ONE of the
  ! tan_th2 values.
  !
  ! m.j. wolff/b.a. whitney, 2000/02/04
  ! history:
  ! 00/02/09 (mjw):  add case to trap theta=90 as indicated by a negative
  ! value of tan_th2.  phiface is called to find the
  ! distance to the theta=90 (xy-plane) surface.  THIS
  ! WILL SLOW DOWN ALL theta<>90 calls as well.
  ! 00/03/16 (mjw):  Using double precision (real(8)).
  ! Define array IND as integer :: (was accidentally real).
  ! Identify only positive roots (t.gt.0)
  ! 09/06/05 (tom):  Catch first and last theta wall cases (don't look for
  ! intersections for 0 and 180 walls)
  !


  implicit none

  ! ... scalar arguments
  real(8) :: ux,uy,uz,x1,y1,z1,tan_th2_1,tan_th2_2

  integer :: i1,i2

  ! ... local scalars
  real(8) :: aa,bb,cc,descr,tmp1,tmp2,tmp3,tan_th2a,tan_th2b

  tmp1 = x1*ux + y1*uy
  tmp2 = ux*ux + uy*uy
  tmp3 = x1*x1 + y1*y1

  !
  ! put minimum value of tan_th2 in 'a' variable, and switch the order
  ! of the "face" offsets if necessary.  This makes it easier to test
  ! for possible occurance of theta=90 (as flagged by -1. value)
  !
  if (tan_th2_2.lt.tan_th2_1) then
     tan_th2a = tan_th2_2
     tan_th2b = tan_th2_1
     ! ... IND array indcates which "face" is associated with which roots
     i1=4
     i2=3
  else
     tan_th2a = tan_th2_1
     tan_th2b = tan_th2_2
     i1=3
     i2=4
  end if

  ! if tan_th2a lt 0, then theta=90 surface is now xy-plane.
  if(tan_th2a.lt.0.d0) then
     ! print*,'trapping theta=90'
     if(uz.ne.0.d0) call insert_t(-z1/uz,i1)
  else
     bb = 2.d0*(tmp1 - tan_th2a*uz*z1)
     aa = tmp2 - tan_th2a*uz*uz
     cc = tmp3 - z1*z1*tan_th2a
     if(abs(aa) > 1.d-10) then
        descr = bb*bb - 4.d0*aa*cc
        if (descr.gt.0.d0.and.tan_th2a.ne.0.d0) then
           descr = sqrt(descr)
           call insert_t((-bb + descr)*0.5d0/aa,i1)
           call insert_t((-bb - descr)*0.5d0/aa,i1)
        end if
     else
        if(abs(bb) > 1.d-10) call insert_t(-cc/bb,i1)
     end if
  end if

  ! ASSUME tan_th2b CAN NEVER be .lt.0
  bb = 2.d0*(tmp1 - tan_th2b*uz*z1)
  aa = tmp2 - tan_th2b*uz*uz
  cc = tmp3 - z1*z1*tan_th2b
  if(abs(aa) > 1.d-10) then
     descr = bb*bb - 4.d0*aa*cc
     if (descr.gt.0.d0.and.tan_th2b.ne.0.d0) then
        descr = sqrt(descr)
        call insert_t((-bb + descr)*0.5d0/aa,i2)
        call insert_t((-bb - descr)*0.5d0/aa,i2)
     end if
  else
     if(abs(bb) > 1.d-10) call insert_t(-cc/bb,i2)
  end if
  return
end subroutine thetaface
