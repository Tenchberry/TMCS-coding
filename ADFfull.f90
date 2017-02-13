!Program to find the angular distribution function of a periodic lattice
!Will take user input and calculate all distances between atoms to find the
!bond angles between all triplets of atoms within the system which are then used
!to form a probability distribution function.
!Last Modified: 3/2/2017
!Khaled Abdel Maksoud

program adffull
  implicit none

  integer, parameter :: maxatoms = 4000
  real, dimension (maxatoms,3) :: R
  real, dimension (0:180) :: bin, final_bin
  real, dimension (3) :: total, a, b, c
  character (len=80) :: input, outfile, atomtype
  real :: magA, magB, magC, magA2, magB2, magC2
  real :: cosine1, cosine2, angle1, angle2
  real :: thickness
  integer :: natoms, ntype, nxcells, nycells, nzcells, nangles
  integer :: i, j, k, l, m, n, u, v
  real :: pi = acos(-1.0)

  bin(0:180) = 0

!Opening files - next step is to enter different subroutines for reading pdb and xyz files

  100 write (6,*) 'Enter output filename for full ADF: '
  read (*,200) outfile
  200 format (a)
  open (unit=200, file=outfile, status='new', err=100)

  300 write (6,*) 'Enter input filename: '
  read (*,400) input
  400 format (a)
  open(unit=400, file=input, status='old', err=300)

!Reading lattice information from inputted datafile

  read(400,*) natoms
  read(400,*) ntype
  read(400,*) nxcells, nycells, nzcells
  read(400,*) total(1), total(2), total(3) !Width variables will hold the total distance covered by the periodic system
  thickness = total(1)/nxcells

!Inputting positions - array R will store the positions r between atoms i and j

   do 50 i = 1,natoms
     read (400,*) atomtype, R(i,1), R(i,2), R(i,3)
   50 continue

!Main angular distribution statements

  nangles = (natoms*(natoms-1)*(natoms-2))/2

  do 60 i = 1,natoms-1
    do 70 j = i+1,natoms

      call distances(magA, magA2, a, j, i)

      do 90 k = 1,natoms
        if (k.ne.i.and.k.ne.j) then

          call distances(magB, magB2, b, j, k)
          call distances(magC, magC2, c, k, i) !All distances now collected to calculate bond angles

          call histogram(magA, magB, a, b, cosine1, angle1, u)
          call histogram(magA, magC, c, a, cosine2, angle2, v)

        end if
      90 continue
    70 continue
  60 continue

!Generating output and closing files

  write(*,*) 'Angular distribution function succesfully calculated, Normalising...'

  do 900 m=0,180
    final_bin(m) = bin(m)/nangles !dividing total bin count by nangles to obtain probablities for each angle between 0 to 180
    write(200,*) m, final_bin(m) !Plots probablities for each angle into the output file
  900 continue

  endfile (unit=200)
  close (unit=200) !Ending and closing outputfile

contains
  subroutine distances (magX, magX2, x, val1, val2)

    !Calculating bond lengths between all atoms

    !Local variables
    real :: magX, magX2
    real, dimension(3) :: x
    integer :: val1, val2

    magX2 = 0.0
    do 80 l = 1,3
      x(l) = R(val1,l) - R(val2,l)       !Calculating bond length
      if (abs(x(l)).gt.(total(l)/2)) then
        if (x(l).gt.0.0) then
          x(l) = x(l) - total(l)         !Applying periodic boundary conditions
        else
          x(l) = x(l) + total(l)
        end if
      end if
      magX2 = magX2 + x(l)**2
    80 continue
    magX = sqrt(magX2)

  end subroutine distances

  subroutine histogram (magY1, magY2, r1, r2, cosx, angle, intangle)

    !Calculating bond angles and binning angles into bins of histogram to generate angular distribution function

    !Local variables
    real :: magY1, magY2, cosx, angle
    integer :: intangle
    real, dimension(3) :: r1, r2

    cosx = (r1(1)*r2(1) + r1(2)*r2(2) + r1(3)*r2(3))/(magY1*magY2) !Bond angle calculated
    if (cosx.gt.1.0) then
      cosx = 1.0
    elseif (cosx.lt.-1.0) then
      cosx = -1.0              !Scales angles such that the maximum angle is 180 degrees, minimum is 0 degrees
    end if
    angle = acos(cosx)*180.0/pi !angle (in degrees) ijk
    intangle = nint(angle)
    bin(intangle) = bin(intangle) + 1 !adds angle to bincounter of angles in ADF

  end subroutine histogram

end program adffull
