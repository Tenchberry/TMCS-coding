!Program to find the angular distribution function of a periodic lattice
!Finds an inner distribution function within a radial cutoff and an
!outer distribution between atoms pairs within the cutoff boundary and their nearest neighbours out of the boundary
!Last Modified: 3/2/2017
!Khaled Abdel Maksoud

program adf
  implicit none

  integer, parameter :: maxatoms = 4000
  real, dimension (maxatoms,3) :: R
  real, dimension (0:180) :: binout, binin
  real, dimension (3) :: total, a, b, c
  character (len=80) :: input, outfile1, outfile2, atomtype
  real :: magA, magB, magC, magA2, magB2, magC2
  real :: weight1, weight2, cosine1, cosine2, angle1, angle2
  real :: shell, shellinit, thickness
  integer :: natoms, ntype, nxcells, nycells, nzcells
  integer :: i, j, k, l, m, n, u, v
  real :: pi = acos(-1.0)

  binout(0:180) = 0
  binin(0:180) = 0  !Initialising histogram bincounts for inner ADF (binin) and outer ADF (binout)

!Opening files - next step is to enter different subroutines for reading pdb and generated xyz files

  55 write(6,*) 'Enter output filename for outer ADF: '
  read (*,75), outfile1
  75 format (a)
  open (unit=75, file=outfile1, status='new', err=55)
  100 write (6,*) 'Enter output filename for inner ADF: '
  read (*,200) outfile2
  200 format (a)
  open (unit=200, file=outfile2, status='new', err=100)

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

  write(6,*) 'Enter radial cutoff distance (angstroms): '

!Recommendation for radial cutoff distance given for the lattice type of the system - based on radial distribution function calculation

  if (ntype.eq.1) then
    write(*,*) 'Simple cubic lattice detected - recommended radial cutoff = 1.2'
  elseif (ntype.eq.2) then
    write(*,*) 'Face-centred cubic lattice detected - recommended radial cutoff = 1.2'
  elseif (ntype.eq.3) then
    write(*,*) 'Body-centred cubic lattice detected - recommended radial cutoff = 1.4'
  elseif (ntype.eq.4) then
    write(*,*) 'Diamond-type lattice detected - recommended radial cutoff = 1.3'
  elseif (ntype.eq.5) then
    write(*,*) 'Ideal gas lattice model detected - recommend setting radial cutoff to larger than the system lattice size'
  else
    write(*,*) 'Lattice type is unidentified'
  end if

  read(5,*) shellinit
  shell = thickness*shellinit !corrects the radial cutoff distance for the size of the unit cell

!Inputting positions - array R will store the positions r between atoms i and j

   do 50 i = 1,natoms
     read (400,*) atomtype, R(i,1), R(i,2), R(i,3)
   50 continue

!Main angular distribution statements

  do 60 i = 1,natoms-1
    do 70 j = i+1,natoms

      call distances(magA, magA2, a, j, i)

      if (magA.le.shell) then             !Distances between atoms within the radial cutoff found
        do 90 k = 1,natoms
          if (k.ne.i.and.k.ne.j) then

            call distances(magB, magB2, b, j, k)
            call distances(magC, magC2, c, k, i)

            call histogram_partition(magA, magB, a, b, weight1, cosine1, angle1, u)
            call histogram_partition(magA, magC, c, a, weight2, cosine2, angle2, v)

          end if
        90 continue
      else
        go to 60
      end if
    70 continue
  60 continue

!Generating output and closing files

  write(*,*) 'Outer ADF - angular distribution function succesfully calculated'
  write(*,*) 'Inner ADF - angular distribution function succesfully calculated'

  do 900 m=0,180
    write(75,*) m, binout(m) !Writing outer ADF histogram bincounts into output
    write(200,*) m, binin(m) !Writing inner ADF histogram bincounts into output
  900 continue

  endfile (unit=200)
  close (unit=200)
  endfile (unit=75)
  close (unit=75)  !Ending and closing outputfiles for inner and outer ADF

contains
  subroutine distances (magX, magX2, x, val1, val2)

    !Calculating bond lengths between all atoms

    !Local variables
    real :: magX, magX2
    real, dimension(3) :: x
    integer :: val1, val2

    magX2 = 0.0
    do 80 l = 1,3
      x(l) = R(val1,l) - R(val2,l)
      if (abs(x(l)).gt.(total(l)/2)) then
        if (x(l).gt.0.0) then
          x(l) = x(l) - total(l)
        else
          x(l) = x(l) + total(l)
        end if
      end if
      magX2 = magX2 + x(l)**2
    80 continue
    magX = sqrt(magX2)

  end subroutine distances

  subroutine histogram_partition (magY1, magY2, r1, r2, weightx, cosx, angle, intangle)

    !Calculating bond angles and binning angles into histogram to generate angular distribution function

    !Local variables
    real :: magY1, magY2, weightx, cosx, angle
    integer :: intangle
    real, dimension(3) :: r1, r2

    weightx = (magY1*magY2)**4 !weighting factor for prioritisation of bond angles between nearest neighbours
    cosx = (r1(1)*r2(1) + r1(2)*r2(2) + r1(3)*r2(3))/(magY1*magY2)
    if (cosx.gt.1.0) then
      cosx = 1.0
    elseif (cosx.lt.-1.0) then
      cosx = -1.0
    end if
    angle = acos(cosx)*180.0/pi !angle (in degrees) ijk
    intangle = nint(angle)
    if (magY2.gt.shell) then
      binout(intangle) = binout(intangle) + 1/weightx !adds angle to bincounter of angles in outer ADF
    else
      binin(intangle) = binin(intangle) + 1/weightx !adds angle to bincounter of angles in inner ADF
    end if

  end subroutine histogram_partition

end program adf
