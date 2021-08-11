                     program color
c$$$$ calls getmat getone major peruse
c
c  Main program for color contouring of rectangular arays.
c  Program runs under a simple interpretive language with commands
c  like
c    file topography
c    read 43 55
c  See color.doc for a complete description.
c  The program generates a PostScript file for plotting
c
c  The principle of the program is to describe all the boundaries of
c  the different colored regions (the contour lines) which are forced
c  to be closed curves.  These closed polygons are filled with the 
c  appropriate color using the PostScript "fill" command.  Of course
c  this must be done in a proper order to prevent covering up previously
c  drawn regions; the program plots in order of decreasing regional
c  area
c
c  Language interpretation is handled by the routines  peruse,
c  getval,  getchr.  Peruse reads in all commands up to certain
c  key points (read or plot), saving the list in a character array
c  for subsequent analysis by getval and getchr.  These routines
c  are called whenever values or strings are needed at the appropriate
c  point in the program
c
c  There are essentially 7 separate program units:
c  (1) the main program, start up and quit routines.
c  (2) the language unit                                    (peruse.f)
c  (3) the array input routines                 (getmat.f)
c  (4) the map generator and page manager       (major.f)
c  (5) the contour line generating algorithm    (glevel.f)
c  (6) the level-line output handler                        (dumpit.f)
c  (7) the color palette data base              (artist.f)
c
      parameter (inmx=200)
      character *80 input(inmx), output
      common /dicta/ input
      common /outfil/ output
      common /ndict/ iecho,nin
c
      parameter (maxa=920 000)
      common // a(maxa)
c
      common /xtreme/ bord,amin,amax,alim(3)
      common /inout/ inp,iout,inf,info
c
      parameter (levmx=100)
      common /noir/ scutum,kolors,art(4,levmx),kart(levmx),loutln(levmx)
c
      common /atlas/ map,mtab,ntab,xtab(201),ytab(201)
c
      data ihead/1/, m,n/0,0/
c
ce
      iecho=0
c  
c  
      do 5000 k=1, inmx
c  Read commands up to 'plot' or 'quit'
        call peruse(ihead)
        call getone('verbose', verb, nf)
        if (nf .eq. 1) info=nint(verb)
c
c  Input the array
        call getmat(m, n, a)
        ihead=0
c
c  Draw the map
        call major(m+2, m, n, a)
        ihead=1
 5000 continue
      end
c______________________________________________________________________
      subroutine quit
c$$$$ calls getone system encaps
c  Brings program to an orderly conclusion;
c  Send file to screen if requested.
c  CAUTION: "system" is a SUN routine to invoke a UNIX shell and
c           show is local window display utility.  Users on other
c           operating systems may delete the call or provide
c           an alternate
      character*80 output
      common /outfil/ output
      common /inout/ inp,iout,inf,info
      data d1,d2,d3,d4/0,0,0,0/
c
      write(9, '(a)')'showpage'
      call encaps(3, d1,d2,d3,d4)
      close (unit=9)
      call getone('show', dum, nshow)
c  CAUTION: Next line is system-dependent.  Delete call to system
c  if not supported in your compiler.  gv is ghostview in disguise
      if (nshow .ge. 0) call system
     $('gv '//output//'&')
      if (output .ne. 'mypost') write(*,'(a,a)')
     $' PostScript file written to: ',output
      write(*,'(a)')' ',' Normal termination'
      stop
      end
c______________________________________________________________________
      blockdata start
c
      common /ndict/ iecho,nin
      common /xtreme/ bord,amin,amax,alim(3)
      common /inout/ inp,iout,inf,info
c
      data iecho/0/, nin/0/
      data alim/0,0,0/
c
c inp=std input; iout=std output; inf=array diskfile; info=print freq
      data inp/5/,       iout/6/,             inf/7/, info/1/
c
      end
c______________________________________________________________________
*====================================
* Level-line output handler   BEGIN===================
*
      subroutine dumpit(val, nxy, xx, yy)
c$$$$ calls mapper
c  Saves in common /arcs/ the vector kxy, which contains closed
c  contour loops as they are generated; also saves path area and
c  inside color value.  Calls mapper to transform saved lines to
c  their new coordinate plane
c
c     Arguments
c  val    is the current contour level
c  nxy    is the number of points defining the contour, which is
c         always closed
c  xx,yy are arrays containing the points on the contour; these are
c        such that a grid step is one unit in each, (0,0) is the top
c        left corner, x increases with i, the row counter, and y with
c        j the column counter
c
      dimension xx(*),yy(*)
      parameter (lpmx=2500, knmx=50 000)
      integer*2 kxy
      common /perim/ ltrue(lpmx)
      common /arcs/ loops,istart(lpmx),area(lpmx),lev(lpmx),kxy(2,knmx)
      common /trace/ level,em,en,xfac,yfac
c
      parameter (maxa=920 000)
      common // a(maxa)
c
      common /inout/ inp,iout,inf,info
c
c  Advance the loop counter; save next start position; trap errors
      loops=loops+1
      if (loops .gt. lpmx) then
        write(*,'(a)') ' ',
     $  '>>> Sorry - too many contour lines.  Array may be transposed.',
     $  'If not, reduce number of levels, or recompile color with',
     $  'increased parameter LPMX'
        stop
      endif
c
c  Calculate -2*area of this loop
      surf=0.0
      xx(nxy+1)=xx(1)
      yy(nxy+1)=yy(1)
      do 1100 i=2, nxy+1
        surf=surf + (xx(i)*yy(i-1) - xx(i-1)*yy(i))
 1100 continue
      area(loops)=-abs(surf)
c
c  Discover which color lies INSIDE contour: midf= +1 if high values
c  inside, -1 for low values.
c  Move along contour until one of coordinates lies on a grid line
      do 1200 i=1, nxy-1
        if (nint(xx(i)).eq.xx(i) .or. nint(yy(i)).eq.yy(i)) goto 1210
 1200 continue
      i=1
 1210 i1=nint(xx(i) + 1.0)
      j1=nint(yy(i) + 1.0)
      ij=(j1-1)*(em+2.0) + i1
      midf=sign(1.0, ((yy(i+1)-yy(i))*(i1-1.0- xx(i)) -
     $(xx(i+1)-xx(i))*(j1-1.0 - yy(i)))*(a(ij) - val)*surf)
      lev(loops)=max(1, level + (midf - 1)/2)
      ltrue(loops)=level
c
c  Convert to thousandths of an inch and store in common /arcs/;
c  constrain path to lie within or on boundary of inner rectangle
      j=istart(loops)
      call mapper(kxy(1,j), xx(1), yy(1))
      j=j+1
      call mapper(kxy(1,j), xx(2), yy(2))
c  From this point on, remove long runs of constant x or y values which
c  arise in the border
      do 1500 i=3, nxy
        j=j+1
        call mapper(kxy(1,j), xx(i), yy(i))
        if (kxy(1,j).eq.kxy(1,j-2) .and. kxy(1,j).eq.kxy(1,j-1)) then
          kxy(1,j-1)=kxy(1,j)
          kxy(2,j-1)=kxy(2,j)
          j=j - 1
          goto 1500
        endif
        if (kxy(2,j).eq.kxy(2,j-2) .and. kxy(2,j).eq.kxy(2,j-1)) then
          kxy(1,j-1)=kxy(1,j)
          kxy(2,j-1)=kxy(2,j)
          j=j - 1
          goto 1500
        endif
 1500 continue
c
      istart(loops+1)=j + 1
      if (istart(loops+1) .ge. knmx) then
        write(*,'(a)') ' ',
     $  '>>> Sorry - total length of contour lines is too long.',
     $  'Array may be transposed.  If not, use fewer contour levels,',
     $  'or recompile color with bigger parameter KNMX'
        stop
      endif
c
      return
      end
c_______________________________________________________________________
      blockdata mappng
      common /atlas/ map,mtab,ntab,      xtab(201),  ytab(201)
      data       map/1/, mtab,ntab/0,0/, xtab/201*1/,ytab/201*1/
      end
c_______________________________________________________________________
      subroutine mapper(kxy, xp, yp)
c$$$$ calls aitoff
c  Maps the point (xp, yp) into 1/1000 inch scaled version.
c  Values of xp are forced into the interval (1,em) and
c  yp into (1,en);  the scale factors and limits common through /trace/
c
c  Choice of map comes through /atlas/ in map
      integer*2 kxy
      dimension kxy(*)
      common /trace/ level,em,en,xfac,yfac
      common /atlas/ map,mtab,ntab,xtab(201),ytab(201)
c
c  Constrain coordinate to lie within the true array, sans border
      x=max(1.0, min(em, xp))
      y=max(1.0, min(en, yp))
c
c  Select the mapping
      goto (1000, 2000, 3000), map
c
c  The identity mapping for our purposes.  Note the conversion to
c  integer*2 thousandths of an inch
 1000 kxy(1)=xfac*x
      kxy(2)=yfac*y
      return
c
c  Aitoff equal-area from presumed scaled long-lat input pairs
 2000 alat =179.9*(y-1.0)/(en-1.0) - 89.95
      along=359.9*(x-1.0)/(em-1.0) -179.95
      call aitoff(alat, along, x, y)
      x=1.0 + (em-1.0)*(0.25*x + 0.5)
      y=1.0 + (en-1.0)*(0.50*y + 0.5)
      goto 1000
c
c  Tensor-producr representation from table in /atlas/
 3000 n=x
      x=xtab(n) + (x-n)*(xtab(n+1) - xtab(n))
      n=y
      y=ytab(n) + (y-n)*(ytab(n+1) - ytab(n))
      goto 1000
c
      end
c_______________________________________________________________________
      subroutine aitoff(lat, lon, x, y)
c$$$$ calls nothing
c  Maps latitude, longitude into Aitoff's equal-area projection.
c  The sphere is mapped into an ellipse, major semi-axis 2, minor
c  semi-axis 1; the equator is imaged into the major axis, the
c  Greenwhich meridian into the minor axis.  The origin (0, 0),
c  where the axes cross, may be moved to longitude  phi0  by
c     call aitoff(lat, lon-phi0, x, y)
c
      real lat, lon
      data degby2/0.008726646/, deg/0.0174532925/
c
      clat=cos(deg*lat)
      half=degby2*(mod(lon+540.0, 360.0) - 180.0)
      R=sqrt(1.0/(abs(1.0 + clat*cos(half)) + 1.0e-15))
      y=R*sin(deg*lat)
      x=2.0*R*clat*sin(half)
      return
      end
c__________________________________________________________________
      subroutine tensor
c$$$$ calls getchr getone
c  Read spacing file and set up table
      common /atlas/ map,mtab,ntab,xtab(201),ytab(201)
      common /trace/ level,em,en,xfac,yfac
      common /inout/ inp,iout,inf,info
      character *80 space
c
c  Choose a mapping function
      call getone('mapping', amap, nap)
      if (nap .eq. 1) map=amap
c
c  See if tensor spacing is specified
      call getchr('tensor', space, ntens)
c  No tensor spacing mentioned
      if (ntens .lt. 0) then
        if (map .eq. 3 .and. (ntab.ne.en .or. mtab.ne.em)) then
          write(*,'(a)')
     $    '>>> Tensor spacing and array size are inconsistent:',
     $    ' Map drawn with even spacing'
          map=1
        endif
c  Read from file if requested.
c  Read list of x spacings first, then y spacings
      elseif (ntens .eq. 0) then
        write(*,'(a)') '>>> You must state a filename in the command'
      elseif (ntens .ge. 1) then
c  Open a new file
        open (unit=8, file=space,err=5000)
        mtab=em
        xtab(1)=0.0
        read (8,*, err=5010) (xtab(j),j=2, mtab)
        ntab=en
        ytab(1)=0.0
        read (8,*, err=5010) (ytab(j),j=2, ntab)
        map=3
        close (unit=8)
c  Sum increment and construct stored table
        do 1100 j=2, mtab
          xtab(j)=xtab(j) + xtab(j-1)
 1100   continue
        do 1150 j=1, mtab
          xtab(j)=(em-1.0)*xtab(j)/xtab(mtab) + 1.0
 1150   continue
c   Do the same for y
        do 1200 j=2, ntab
          ytab(j)=ytab(j) + ytab(j-1)
 1200   continue
        do 1250 j=1, ntab
          ytab(j)=(en-1.0)*ytab(j)/ytab(ntab) + 1.0
 1250   continue
      endif
      if (map .ne. 1) write(*,'(/a,i3)')
     $' Color map is drawn under mapping ',map
      return
c
c  Error returns
 5000 write(*,'(a,a)') '>>> Unable to open file ',space
      return
c
 5010 write(*,'(a,a)') '>>> Read error in file ', space
      return
      end
c__________________________________________________________________
* =================== Level-line output handler   END
*====================================
*====================================
* Contour line generating algorithm BEGIN================
*
      subroutine glevel(mdim, m, n, f, ds, zlevel)
c$$$$ calls circ dotted draw1 draw2 halt parab
c
c  Principal routine for contouring rectangular arrays of real data.
c  Algorithm by  Dr David P. Anderson,  modified by R. L. Parker (1984).
c
c  Arguments
c   mdim   1st dimension of  f in caller
c   m      number of rows of data in array  f
c   n      number of columns
c   f      2-dimensional array of values to be contoured
c   ds     step between points on contour in cell units
c   zlevel contour level traced out in this call
c
c  Conventions
c  The standard plotting convention is that  x  increases with
c  i,  the row index of  f(i,j)  and  y  increases with   j.
c  The point (0, 0) with respect to current plot origin
c  corresponds to the matrix element  f(1,1).
c
c  Notes
c  The system calls  dumpit  which has the responsibility of
c  disposing of a vector of points along each contour.  This may
c  mean plotting the points or saving them in an array for later use.
c
c  common /inout/ controls input and output.
c  inp,iout  the input and output unit numbers of the system.
c            Actually glevel nevers uses  inp
c
      parameter (maxpts=998, maxcmp=10*maxpts)
      logical*1 flag
      dimension f(mdim,*)
      common /hlevel/ value
      common /inout/ inp,iout,inf,info
      common /xtreme/ bord,amin,amax,alim(3)
      common /option/ smooth,
     $                xvals(maxpts),yvals(maxpts)
      common /big1/ point(2,maxpts,maxpts) /big2/ flag(2,maxpts,maxpts)
      common /little/ iline,eps,eps2,bound(4),x(maxcmp),y(maxcmp)
c
      integer off(2,3,2),dir
      real c1(2),c2(2),bx(3),by(3),p1(2),p2(2)
      real ang1(2),ang2(2),ang3(2),ang4(2),ang5(2),ang6(2)
      logical ifcase(3),closed
      integer xoff(3,2,2),yoff(3,2,2),koff(3,2,2),doff(3,2,2)
      data xoff/0,0,1,0,1,0, 0, 0, 1,-1,-1,-1/
      data yoff/0,1,0,1,0,0,-1,-1,-1, 0, 0, 1/
      data koff/2,1,2,1,2,1, 2, 1, 2, 1, 2, 1/
      data doff/2,1,1,1,1,2, 2, 2, 1, 2, 2, 1/
      data off/-1,0,1,0,2,0,0,-1,0,1,0,2/
      data ix0, iy0, ik0/0,0,0/
c
c  Catch various errors
c
      if (m.gt.maxpts .or. n.gt.maxpts) call
     $ halt ('Array dimensions too large - increase parameter  MAXPTS')
c
c  Load  xvals, yvals  with coordinates in inches if dx, dy not zero
c
      do 231 i=1, m
        xvals(i)=i-1
 231  continue
      do 232 j=1, n
        yvals(j)=j-1
 232  continue
      smooth=1.0/ds
c
c
      bound(1)=xvals(1)
      bound(2)=xvals(m)
      bound(3)=yvals(1)
      bound(4)=yvals(n)
      eps2=1e-5*max(bound(4)-bound(3), bound(2)-bound(1))**2
c
c
c
c  Largest, smallest function values from /xtreme/
c
      range=amax - amin
c  Trace the contour level
c
        val=zlevel
        value=val
c
c  Make sure no grid points have exact contour values
      if (range .eq. 0.0) return
      tol=1e-7
      do 1005 j=1, n
        do 1006 i=1, m
            if (abs(f(i,j)-zlevel) .le. abs(zlevel)*tol) then
              if (zlevel.ne. 0.0) f(i,j)=zlevel*(1.0+ tol)
              if (zlevel.eq. 0.0) f(i,j)=range*tol
            endif
 1006   continue
 1005 continue
c
c
c  Find the points where contour intersects grid lines
c
        do 1009 j=1, n
          do 1010 i=1, m
            a=f(i,j)
            do 1011 k=1, 2
              i2=i+off(1,2,k)
              j2=j+off(2,2,k)
c
c  Make sure segment isn't off the grid
c
              if (i2.le.m.and.j2.le.n) then
                b=f(i2,j2)
                if (a.lt.val.and.b.gt.val.or.b.lt.val.and.a.gt.val)then
c
c  Here if contour crosses this segment - set flag and find point
c
                  flag(k,i,j)=.true.
                  i1=i+off(1,1,k)
                  j1=j+off(2,1,k)
                  i3=i+off(1,3,k)
                  j3=j+off(2,3,k)
c
c  If segment is on edge of grid (pointing toward edge)
c  use single-parabola approximation
c
                  if (i1.le.0.or.j1.le.0) then
                    c=1.-parab(f(i3,j3),f(i2,j2),f(i,j),val)
                    point(k,i,j)=c
                  else if (i3.gt.m.or.j3.gt.n) then
                    c=parab(f(i1,j1),f(i,j),f(i2,j2),val)
                    point(k,i,j)=c
                  else
c
c  Otherwise use double-parabola approximation
c
                    c=parab(f(i1,j1),f(i,j),f(i2,j2),val)
                    d=1.-parab(f(i3,j3),f(i2,j2),f(i,j),val)
                    point(k,i,j)=(c+d)/2.0
                  end if
                else
                  flag(k,i,j)=.false.
                end if
              else
                flag(k,i,j)=.false.
              end if
1011        continue
1010      continue
1009    continue
c
c  All points found -- now look for components
c
1012    continue
c
c  The variables in this part are
c  ix, iy           indices of start of current edge
c  ik               direction of edge, 1=x, 2=y
c  dir              direction out of edge, 1=up/right, 2=left/down
c  closed           true if component is a loop
c  ix0,iy0          if closed component, starting point
c
c  Search for a point on border to start from
c
          ik=1
          closed=.false.
          do 1013 ix=1, m-1
            if (flag(1,ix,1)) then
              iy=1
              dir=1
              goto 100
            else if (flag(1,ix,n)) then
              iy=n
              dir=2
              goto 100
            end if
1013      continue
          ik=2
          do 1014 iy=1, n-1
            if (flag(2,1,iy)) then
              ix=1
              dir=1
              goto 100
            else if (flag(2,m,iy)) then
              ix=m
              dir=2
              goto 100
            end if
1014      continue
c
c  If none, search for point in middle
c
          do 1015 iy=1, n
            do 1016 ix=1, m
              do 1017 ik=1, 2
                if (flag(ik,ix,iy)) then
                  dir=1
                  closed=.true.
                  ix0=ix
                  iy0=iy
                  ik0=ik
                  goto 100
                end if
1017          continue
1016        continue
1015      continue
c
c  If can't find a starting point, we're done with this contour value
c
          goto 1018
c
c  Here when found starting point of a component
c  compute its xy coordinates
c
100       continue
          ax=xvals(ix)
          ay=yvals(iy)
          if (ik.eq.1) then
            ax=ax+point(1,ix,iy)*(xvals(ix+1)-xvals(ix))
          else
            ay=ay+point(2,ix,iy)*(yvals(iy+1)-yvals(iy))
          end if
c
c  Start a list (in x, y) of the points on the component
c
          x(1)=ax
          y(1)=ay
          npts=1
          if (.not.closed)flag(ik,ix,iy)=.false.
c
c  Loop for points in this component
c
1019      continue
c
c  From the current edge, there are three possible edges to go to.
c  See which of these edges are on the grid and contain a point
c
            do 1021 icase=1, 3
              ifcase(icase)=.false.
              jx=ix+xoff(icase,ik,dir)
              jy=iy+yoff(icase,ik,dir)
              jk=koff(icase,ik,dir)
              if (jx.gt.0.and.jy.gt.0) then
c
c  Note flag is false for nonexistent edges
c
                if (flag(jk,jx,jy)) then
                  ifcase(icase)=.true.
                  bx(icase)=xvals(jx)
                  by(icase)=yvals(jy)
                  if (jk.eq.1) then
                    bx(icase)=bx(icase)+point(1,jx,jy)
     +                              *(xvals(jx+1)-xvals(jx))
                  else
                    by(icase)=by(icase)+point(2,jx,jy)
     +                              *(yvals(jy+1)-yvals(jy))
                  end if
                end if
              end if
1021        continue
c
c  See if there wasn't a point where one was expected
c  If this happens there is a bug in the program.  Report it to d.p.a.
c
          if (.not.(ifcase(1).or.ifcase(2).or.ifcase(3))) then
            write(*,*)'>>> Glevel returns prematurely! Level=',zlevel
            return
          endif
c
c  If any choice, go across the square
c
            if (ifcase(1).and.ifcase(3)) then
              dist1=(ax-bx(1))**2 + (ay-by(1))**2
              dist3=(ax-bx(3))**2 + (ay-by(3))**2
              if (dist1.lt.dist3) then
                ibest=1
              else
                ibest=3
              end if
            else if (ifcase(1)) then
              ibest=1
            else if (ifcase(3)) then
              ibest=3
            else
              ibest=2
            end if
c
c  Update the coords, etc.
c
            ix=ix+xoff(ibest,ik,dir)
            iy=iy+yoff(ibest,ik,dir)
            dir=doff(ibest,ik,dir)
            ik=koff(ibest,ik,dir)
c
c  Add point to the list
c
            ax=bx(ibest)
            ay=by(ibest)
            npts=npts+1
c
c  Check for array overflow
c
            if (npts.gt.maxcmp)
     $      call halt('Out of room - raise parameter  maxcmp')
            x(npts)=ax
            y(npts)=ay
c
c  Remove the point from grid so you don't use it again
c
            flag(ik,ix,iy)=.false.
c
c  Check for termination of component
c
            if (closed) then
              if (ix.eq.ix0.and.iy.eq.iy0.and.ik.eq.ik0) then
                npts=npts-1
                goto 1020
              end if
            else
              if (ix.eq.1.and.ik.eq.2.or.
     +        ix.eq.m.and.ik.eq.2.or.
     +        iy.eq.1.and.ik.eq.1.or.
     +        iy.eq.n.and.ik.eq.1) then
                goto 1020
              end if
            end if
            goto 1019
1020      continue
c
c  Here when finished finding points in component.
c  Points are in arrays x,y.  "closed" tells component type
c
          if (.not.closed) then
c
c  Here to draw open curve
c
            if (npts.eq.2) then
c
c  Component has only two points -- just draw line segment
c
              p1(1)=x(1)
              p1(2)=y(1)
              p2(1)=x(2)
              p2(2)=y(2)
              call dotted(p1,p2)
            else
c
c  Component has three or more points -- draw nice curves
c
              call circ(1,2,3,c1,rad1,ang1,ang2,ang3)
              call draw1(c1,rad1,ang1,ang2)
              do 1030 i=2, npts-2
                call circ(i,i+1,i+2,c2,rad2,ang4,ang5,ang6)
                call draw2(c1,c2,rad1,rad2,ang2,ang3,ang4,ang5)
                c1(1)=c2(1)
                c1(2)=c2(2)
                rad1=rad2
                ang2(1)=ang5(1)
                ang2(2)=ang5(2)
                ang3(1)=ang6(1)
                ang3(2)=ang6(2)
1030          continue
              call draw1(c1,rad1,ang2,ang3)
            end if
          else
c
c  Here to draw closed curve
c
            call circ(npts,1,2,c1,rad1,ang1,ang2,ang3)
            do 1031 i=1, npts
              j=mod(i,npts)+1
              k=mod(i+1,npts)+1
              call circ(i,j,k,c2,rad2,ang4,ang5,ang6)
              call draw2(c1,c2,rad1,rad2,ang2,ang3,ang4,ang5)
              c1(1)=c2(1)
              c1(2)=c2(2)
              rad1=rad2
              ang2(1)=ang5(1)
              ang2(2)=ang5(2)
              ang3(1)=ang6(1)
              ang3(2)=ang6(2)
1031        continue
          end if
          goto 1012
1018    continue
c
c  Flush the plot buffer.  The coordinates given cannot be on a contour
c
      iline=0
      c1(1)=xvals(1) - 0.1
      c1(2)=c1(1)
      call dotted(c1, c1)
      return
      end
c______________________________________________________________
      subroutine circ(i, j, k, center, radius, ang1, ang2, ang3)
c$$$$ calls nothing
      parameter (maxpts=998, maxcmp=10*maxpts)
      logical*1 flag
      common /inout/ inp,iout,inf,info
      common /option/ smooth,
     $                xvals(maxpts),yvals(maxpts)
      common /big1/ point(2,maxpts,maxpts) /big2/ flag(2,maxpts,maxpts)
      common /little/ iline,eps,eps2,bound(4),x(maxcmp),y(maxcmp)
c
c  Routine to find the center, radius and angular limits of the
c  circular arc determined by the points with indices i, j, k in the
c  arrays x and y.  If colinear return radius=-1 and put the
c  three points in ang1, ang2 and ang3
c
      real center(2),ang1(2),ang2(2),ang3(2)
      data pi2/6.283185/
c
c  The formulas here are hard to derive, at least they were for me
c
      a=(x(i)+x(j))/2.
      b=(y(i)+y(j))/2.
      c=y(j)-y(i)
      d=x(i)-x(j)
      e=(x(j)+x(k))/2.
      ff=(y(j)+y(k))/2.
      g=y(k)-y(j)
      h=x(j)-x(k)
      disc=d*g-c*h
c
c  Check for colinearity
c
      if (abs(disc).lt.eps2) then
        radius=-1.
        ang1(1)=x(i)
        ang1(2)=y(i)
        ang2(1)=x(j)
        ang2(2)=y(j)
        ang3(1)=x(k)
        ang3(2)=y(k)
      else
        t=(d*a+ff*c-b*c-d*e)/disc
        center(1)=e+g*t
        center(2)=ff+h*t
        radius=sqrt((center(1)-x(i))**2+(center(2)-y(i))**2)
        a=x(i)-center(1)
        b=y(i)-center(2)
        c=x(j)-center(1)
        d=y(j)-center(2)
        e=x(k)-center(1)
        ff=y(k)-center(2)
        ang1(1)=atan2(b,a)
        ang2(1)=atan2(d,c)
        ang3(1)=atan2(ff,e)
c
c  Make sure angles are in the right order by adding pi*2 to
c  some of them, if needed
c
        if (ang1(1).lt.ang2(1)) then
          if (ang1(1).lt.ang3(1).and.ang3(1).lt.ang2(1)) then
            ang1(1)=ang1(1)+pi2
          else if (ang3(1).lt.ang1(1)) then
            ang3(1)=ang3(1)+pi2
          end if
        else
          if (ang2(1).lt.ang3(1).and.ang3(1).lt.ang1(1)) then
            ang1(1)=ang1(1)-pi2
          else if (ang3(1).gt.ang1(1)) then
            ang3(1)=ang3(1)-pi2
          end if
        end if
      end if
      return
      end
c______________________________________________________________
      subroutine draw1(center, radius, ang1, ang2)
c$$$$ calls dotted vmov2
c  Routine to draw a circular arc
c  radius = 0 means draw segments from point ang1 to ang2
c
      parameter (maxpts=998, maxcmp=10*maxpts)
      logical*1 flag
      common /inout/ inp,iout,inf,info
      common /option/ smooth,
     $                xvals(maxpts),yvals(maxpts)
      common /big1/ point(2,maxpts,maxpts) /big2/ flag(2,maxpts,maxpts)
      common /little/ iline,eps,eps2,bound(4),x(maxcmp),y(maxcmp)
      real center(2),radius,ang1(2),ang2(2),p1(2),p2(2)
      data p1/0,0/
c
      if (radius.ge.0) then
        dang=ang2(1)-ang1(1)
        arclng=abs(dang)*radius
        nsteps=arclng*smooth+1
        da=dang/nsteps
        do 1000 i=0, nsteps
          angle=ang1(1)+i*da
          p2(1)=radius*cos(angle)+center(1)
          p2(2)=radius*sin(angle)+center(2)
          if (i.ne.0) call dotted(p1,p2)
          call vmov2(p2,p1)
1000    continue
      else
        call dotted(ang1,ang2)
      end if
      return
      end
c______________________________________________________________
      subroutine draw2(cent1, cent2, rad1, rad2, ang1, ang2, ang3, ang4)
c$$$$ calls dotted vmov2
c  Routine to draw the average of two circular arcs
c
c  cent1 and cent2 are the centers
c  rad1 and rad1 are the radii
c  If either radius is zero it means the arc is really a line
c  segment whose endpoints are in (ang1,ang2) or (ang3,ang4)
c
      parameter (maxpts=998, maxcmp=10*maxpts)
      logical*1 flag
      common /inout/ inp,iout,inf,info
      common /option/ smooth,
     $                xvals(maxpts),yvals(maxpts)
      common /big1/ point(2,maxpts,maxpts) /big2/ flag(2,maxpts,maxpts)
      common /little/ iline,eps,eps2,bound(4),x(maxcmp),y(maxcmp)
      real cent1(2),cent2(2),ang1(2),ang2(2),ang3(2),ang4(2)
      real p1(2),p2(2)
      data p1/0,0/
c
c  Find angular ranges (or vector)
c
      dang11=ang2(1)-ang1(1)
      dang12=ang2(2)-ang1(2)
      dang21=ang4(1)-ang3(1)
      dang22=ang4(2)-ang3(2)
c
c  Find approximate length of arc or segment
c
      if (rad1.ge.0) then
        arclng=abs(dang11)*rad1
      else
        arclng=sqrt(dang11**2+dang12**2)
      end if
c
c  Compute number of steps based on that
c
      nsteps=arclng*smooth+1
c
c  Loop for points along curve
c
      do 1000 i=0, nsteps
        t=real(i)/nsteps
c
c  Compute point on first arc
c
        if (rad1.ge.0) then
          angle=ang1(1)+t*dang11
          r1=cent1(1)+rad1*cos(angle)
          r2=cent1(2)+rad1*sin(angle)
        else
          r1=ang1(1)+t*dang11
          r2=ang1(2)+t*dang12
        end if
c
c  Compute point on second arc
c
        if (rad2.ge.0) then
          angle=ang3(1)+t*dang21
          q1=cent2(1)+rad2*cos(angle)
          q2=cent2(2)+rad2*sin(angle)
        else
          q1=ang3(1)+t*dang21
          q2=ang3(2)+t*dang22
        end if
c
c  Take weighted average.  Note weighting coefficients
c
        p2(1)=r1*(1.-t)+q1*t
        p2(2)=r2*(1.-t)+q2*t
        if (i.ne.0) call dotted(p1,p2)
        call vmov2(p2,p1)
1000  continue
      return
      end
c______________________________________________________________
      real function parab(r1, s1, t1, val) 
c$$$$ calls nothing
c  Suppose a quadratic function f has f(-1)=r1, f(0)=s1, and f(1)=t1,
c  with val between s and t.  Where between 0 and 1 does the parabola
c  assume the value "val"?  Rewritten by Bob Parker, 2013.
c
      r=r1-val
      c=s1-val
      t=t1-val
      a=(r+t)/2.0-c
      if (a.eq. 0.0) then 
        parab=c/(c-t)
        return
      endif
      b=(t-r)/2.0
      disc=sqrt(b**2 - 4.0*a*c)
      parab = -(b + sign(disc, b))/(2.0*a)
      if (parab.ge. 0.0 .and. parab.le. 1.0) return
      parab= c/(a*parab)
      return
      end
c______________________________________________________________
      subroutine vmov2(a, b)
c$$$$ calls nothing
      real a(2),b(2)
c
      b(1)=a(1)
      b(2)=a(2)
      return
      end
c______________________________________________________________
      subroutine mabel(val, len, label)
c$$$$ calls nothing
c  For a real number  val  assumed to be a decimal fraction rounded to
c  6 significant figures (or fewer), makes an equivalent character
c  string  label  of minimal length  len
      character*8 form,local*20,label*(*)
      double precision ten, v
      data ten/10.0d00/
c
c  Find the number of significant figures in val
c
      nexp=0
      if (val .ne. 0.0) nexp=int(500.01 + log10(abs(val))) - 500
      v=int(abs(val)/ten**(nexp-5) + 0.5)
      nres=nexp - 5
      do 1100 ipow=0, 4
        if (mod(int(v/ten**ipow+ 0.5), 10) .ne. 0) goto 1110
        nres=nres + 1
 1100 continue
 1110 nfld=3 + max(nexp, -nres, nexp-nres)
c
c  Build a format for writing the string
c
      ndpt=max(0, -nres)
      write(form, '( 2h(f, i2, 1h., i1, 1h))') nfld,ndpt
      if (ndpt .eq. 0) nfld=nfld - 1
c  If format is screwed up try a g format:
      if (index(form,'*') .ne. 0) form='(g12.4)'
c
c  Construct a character string of minimal length for the label
c
      write(local, form) val
c  Remove leading blanks
      is=1
      do 1200 i=1, nfld
        if (local(i:i) .eq. ' ') is=1 + is
 1200 continue
      label=local(is:nfld)
      len=nfld - is + 1
      return
      end
c______________________________________________________________
      subroutine dotted(p1, p2)
c$$$$ calls dumpit
c  Collects a buffer of points to be plotted in the arrays  xx, yy.
c  If a discontinuity is encountered or the buffer is filled, flushes
c  the buffer with routine  dumpit
c
      parameter (maxpts=998, maxcmp=10*maxpts)
      logical*1 flag
      common /hlevel/ value
      common /inout/ inp,iout,inf,info
      common /option/ smooth,
     $                xvals(maxpts),yvals(maxpts)
      common /big1/ point(2,maxpts,maxpts) /big2/ flag(2,maxpts,maxpts)
      common /little/ iline,eps,eps2,bound(4),x(maxcmp),y(maxcmp)
      real p1(2),p2(2),xx(maxcmp),yy(maxcmp)
      save xx,yy,nxy
      data nxy/0/
c
c  Advance the buffer count
      nxy=1 + nxy
c
c  Initialize a new buffer
 1000 if (nxy .eq. 1) then
       xx(nxy)=p1(1)
       yy(nxy)=p1(2)
       nxy=2
       xx(nxy)=p2(1)
       yy(nxy)=p2(2)
       return
      endif
c
c  Add another point to the buffer.  Only done if this line piece
c  is contiguous with the previous one and if the buffer isn't full
      if ((xx(nxy-1)-p1(1))**2+(yy(nxy-1)-p1(2))**2.lt.1e-6
     $   .and. nxy.le.maxcmp)  then
        xx(nxy)=p2(1)
        yy(nxy)=p2(2)
        return
      endif
c
c  Flush the buffer.  Noncontiguous line piece has been encountered
      nxy=nxy - 1
      call dumpit(value, nxy, xx, yy)
c  Re-initialize the buffer
      nxy=1
      if (iline .gt. 0) goto 1000
c
c  If  nline is not positive reinitialize buffer
c
      nxy=0
      return
      end
c______________________________________________________________
      subroutine halt(messag)
c$$$$ calls nothing
c  Program prints the message  messag  then stops
c     common /inout/ inp,iout,inf,info
      character*(*) messag
c
      write(*,'(/1x,2a)') '>>> ',messag
      stop
      end
c______________________________________________________________
*     blockdata dfault
c  Sets various parameters into common for the user's convenience.
c  Values for standard  input and output defined in graphics BLOCKDATA
c     common /inout/ inp,iout,inf,info
c
*     data inp,iout,inf,info /5,6,7,7/
c
*     end
c______________________________________________________________
* ================ Contour line generating algorithm END
*====================================
*====================================
* Language unit BEGIN==========================
*
      subroutine peruse(ihead)
c$$$$ calls quit
c  The command input routine.
c  Reads command lines from the standard input until eof
c  or the command 'plot' or 'read' is encountered.
c  Saves the lines in the character array input for later
c  parsing by getchr and getval.
c  Prints a list of commands for this group of routines
      parameter (inmx=200)
      character *80 input(inmx),line
      common /dicta/ input
      common /ndict/ iecho,nin
      common /inout/ inp,iout,inf,info
c
c
      if (ihead .ne. 0) write(*,'(a)') ' ',
     $' Enter commands for data input and contouring '
c
      do 1500 l=nin+1, inmx
        read(*,'(80a)', end=2500) line
        if (line(1:1) .ne. ' ') then
          nin=nin + 1
          input(nin)=line
          call icheck(line)
        endif
        if (line(1:4) .eq. 'hist') then
          nin=nin - 1
          write(*,'(/a/(i3,a,a75))') 'Command history to date',
     $  (k,': ',input(k)(1:74), k=1, nin)
        endif
        if (line(1:4) .eq. 'plot') goto 2000
        if (line(1:4).eq.'stop' .or. line(1:4).eq.'quit') call quit
 1500 continue
      write(*,'(a,i6,a)')'>>> Command stack full: only',inmx,
     $' allowed','Recompile color with larger parameter inmx'
      stop
c
 2000 return
c
 2500 call quit
      end
c______________________________________________________________________
      subroutine icheck(line)
c$$$$ calls nothing
c  Checks the command fragment com against the catalog; prints a
c  warning if com is not present in the list.
      character*80 line, com*4
c
      com=line(1:4)
      if (0 .eq.
     $+index('abov affi axes dash file form heav heig inte ', com)
     $+index('labe land leve line mapp noba note orie outl ', com)
     $+index('outp pale pall plot read show size skip symb ', com)
     $+index('bord tabl tens verb weig widt quit stop ',com))
     $ write(*,'(a,a)')'>>> Unrecognized command: ',com
      return
      end
c______________________________________________________________________
      subroutine getval(code, values, nwant, nfound)
c$$$$ calls nothing
c
c  Seek most recent instance of the command  code  in the input store.
c  If  code  is not in the store, set nfound=-99.
c  If the command is present, nfound is set  nfound >= 0.
c  Then the remainder of the line is analysed, looking for an array
c  of up to  nwant  numbers.   On return the numbers are put into  the
c  array  values, and  nfound  is set to the number of items found which
c  will not be more than  nwant  but may be zero.
c  If an error is discovered  nfound=-n,  where  n  is the number of
c  numbers properly decoded
c
      parameter (inmx=200)
      character *80 input(inmx),line
      common /dicta/ input
      common /ndict/ iecho,nin
      common /inout/ inp,iout,inf,info
      dimension values(*)
      character local*40, code*4, char*80
c
ce
      do 1010 lin=nin, 1, -1
        line=input(lin)
        if (code .eq. line(1:4)) then
          if (iecho .ge. 1) write(*,'(2a)')'>>> ',line
          nb=index(line, ' ')+1
          char=line(nb:80)
          kn=80 - nb + 1
          goto 1020
        endif
 1010 continue
c  code word not found
      nfound=-99
      return
c
 1020 continue
      k1=1
c  Up to  nwant  numbers are sought
      do 1800 nval=1, nwant
        do 1100 k=k1, kn
          if (char(k:k) .ne. ' ') goto 1200
 1100   continue
        nfound=nval-1
        return
 1200   do 1300 l=k, kn
          if (char(l:l).eq. ',' .or. char(l:l).eq. ' ') goto 1500
 1300   continue
c  Augment unpredictable error trapping of some compilers e.g. Sun
 1500   if (index('     Ee+-.0123456789', char(k:k)) .eq. 0) goto 2000
        k2=l+1
        local=char(k:l-1)
        read (local, '(f40.0)', err=2000) values(nval)
 1800 k1=k2
      nval=1 - nwant
 1900 nfound=1 - nval
      return
c
 2000 write(*, '(a)')' ',
     $'***** Unreadable numbers in this input line:',line
      goto 1900
      end
c______________________________________________________________
      subroutine getchr(code, char, nfound)
c$$$$ calls nothing
c  Seeks a fresh instance of the command  code  in the input store.
c  Previously analysed instances are skipped.
c  If  code  is not present returns  nfound=-99.
c  If the command is present returns the rest of the line, beginning at
c  the first nonblank in  char  and sets nfound=1, unless the line is
c  entirely blank; then  nfound=0  and  char is left undisturbed
c
      parameter (inmx=200)
      character *80 input(inmx),line,code*4, char*(*)
      common /dicta/ input
      common /ndict/ iecho,nin
      common /inout/ inp,iout,inf,info
ce
      do 1010 lin=nin, 1, -1
        line=input(lin)
        if (code .eq. line(1:4)) then
          if (iecho .ge. 1) write(*,'(2a)')'>>> ',line
          goto 1020
        endif
 1010 continue
c  code word not found
      nfound=-99
      return
c
 1020 continue
      nb=index(line, ' ')+1
      do 1200 k=nb, 80
        if (line(k:k) .ne. ' ') then
          char=line(k:80)
          nfound=1
          return
        endif
 1200 continue
c
c  Blank field after code word
      nfound=0
      return
      end
c______________________________________________________________
      subroutine getone(code, value, nfound)
c$$$$ calls getval
c
c  Evaluates a single number from the input store /dicta/ associated
c  with the keyword code.  The answer is returned in  value.
c  nfound  may be 1 if a number is read, 0 if the field after
c  code is blank, -99 if code is absent
c
      dimension values(1)
      character code*4
c
      call getval(code, values, 1, nfound)
      if (nfound .eq. 1) value=values(1)
      return
      end
c______________________________________________________________
      subroutine clear(code)
c$$$$ calls nothing
c  Deletes most recent occurrence of the code in the stack
c
      parameter (inmx=200)
      character *80 input(inmx),line,code*4
      common /dicta/ input
      common /ndict/ iecho,nin
      common /inout/ inp,iout,inf,info
ce
      do 1010 lin=nin, 1, -1
        line=input(lin)
        if (code .eq. line(1:4)) then
          input(lin) = '****'
          return
        endif
 1010 continue
      return
      end
c______________________________________________________________
*=========================== Language unit END
c====================================
*====================================
* Array input routines BEGIN======================
*
      subroutine getmat(m, n, a)
c$$$$ calls getchr getone getval matrab matraf group
c  Routine for reading matrices and re-orienting them on the page.
c  When orient command is invoked defines plotted array as follows:
c  N normal - as in a numerical listing; T transposed array;
c  R 90-deg clockwise turn; RR 2 90-deg turns, etc;
c  V  reflect array in veryical line; H reflect in horizontal line.
c
c  ms, ns  dimensions of array as it was written.
c  klock   number of clockwise quarter-turns.
c  m, n    array dimensions of matrix as actually used in program
c
      parameter (maxa = 920 000)
      common /inout/ inp,iout,inf,info
      character*75  format,file,ops,xfile
      dimension a(*), val(37)
      save format,file,skip,kode1,kode2,acoef,bcoef
      data format(1:4)/'    '/, file/' '/, skip/0.0/
      data kode1/0/, xfile(1:1)/';'/
      data acoef,bcoef/1,0/
c
c  Look for orientation command in list
      call getchr('orie', ops, iorien)
      if (iorien.gt. 0) then
        call group(ops, kode1,kode2)
      endif
c
c  Look for the read statement in the command list.
c  Sept 14,2006: Default orientation changed from val(3)=0 to:
      val(3)=1
      call getval('read', val, 3, nfound)
      call clear('read')
      if (nfound .eq. -99) then
        if (m .gt. 0) then
          write(*, '(a)')'Contouring previously read matrix'
          return
        endif
        write(*, '(a)')' ',
     $  '>>> No read statement in list. Provide one'
        m=0
        return
      elseif (nfound .le. 1) then
        write(*, '(a)')' ',
     $  '>>> read must be followed by array dimensions'
        m=0
        return
      else
        ms=nint(val(1))
        ns=nint(val(2))
        klock=nint(val(3))
c  Revise read specs if orient command has been used
        if (iorien.gt. 0) then
          ms=sign(ms, kode1)
          ns=abs(ns)
          klock=kode2
        endif
      endif
c
c  Impossible array dimensions trapped
      if (abs(ms).le.1.or.abs(ns).le.1.or.abs(ms*ns).gt.maxa) then
        write(*, '(/a,2i5/a)')
     $  '>>> Improper array dimensions',ms,ns,'Array cannot be accepted'
        if (abs(ms*ns).gt.maxa) write(*,'(a,i8)')
     $  '>>> Total array size exceeds maximum permitted, ie ',maxa
        m=0
        return
      endif
c
c  Ascertain file name from file command
      call getchr('file', file, nfl)
      call clear('file')
      if (file(1:1) .eq. ' ') then
        write(*,'(a)')' ',
     $  '>>> No file name has been provided for data'
        m=0
        return
      endif
c  Ascertain format if specified
      call getchr('format', format, ignor)
c
c  Open the file
      if (file .eq. '*') then
        inf=inp
        write(*,'(a)')
     $ '>>> Sorry, file * feature has been removed temporarily'
        m=0
      else
        inf=8
        if (nfl .eq. 0) rewind(unit=8)
        if (format(1:1) .eq. ' ') then
          if (xfile .ne. file) then
            if (inp.ne.inf) close(unit=8)
            open (unit=8, file=file, status='OLD', err=2000)
          endif
          call getone('skip', skip, ignor)
          nskip=nint(skip)
          do 1100 j=1, nskip
            read(8,'(1x)')
 1100     continue
        endif
        if (format(1:1) .eq. 'b' .and. xfile.ne.file)
     $  open (unit=8, file=file,
     $     status='OLD', err=2000, form='UNFORMATTED')
      endif
      xfile=file
c
c  Determine the array dimensions as they really are
c
      if (mod(klock,2) .eq. 0) then
        m=iabs(ms)
        n=iabs(ns)
      else
        m=iabs(ns)
        n=iabs(ms)
      endif
c
c  Code the sign configuration into a single integer
c
      jsign=(5 - isign(2,ns) - isign(1,ms))/2
c
c  Now read the thing in,  either formatted or binary.
c  Notice room for border of 1 cell all round
c
      if (format(1:1).ne.'b') call matraf(m+2,m,n, a(m+4), klock, jsign)
      if (format(1:1).eq.'b') call matrab(m+2,m,n, a(m+4), klock, jsign)
*     if (inf .ne. inp) close (unit=8)
c  Perform affine scaling
      call getval('affine', val, 2, nfound)
      if (nfound .eq. -99) return
      if (nfound .eq. 0)  then
        acoef=1.0
        bcoef=0.0
      elseif (nfound .eq. 2) then
        acoef=val(1)
        bcoef=val(2)
      else
        write(*,'(a)')'>>> Error reading affine - values unchanged'
      endif
      if (acoef.eq.1.0 .and. bcoef.eq.0.0) return
      do 1500 i=2, m+1
        do 1400 j=2, n+1
          ij=(j-1)*(m+2) + i
          a(ij)=acoef * a(ij) + bcoef
 1400   continue
 1500 continue
c
      return
c
c  Can't open the file for some reason
 2000 write(*,'(/a)')'>>> Nonexistent file or read permission denied'
      write(*,'(a)') '>>> Filename: '//file
      m=0
      return
      end
c______________________________________________________________
      subroutine matrab(mdim, m, n, a, klack, jsign)
c$$$$ calls nothing
c  Routine actually to read binary diskfile.
c  It is assumed the file was written by rows with a statement like
c        do 1000 i=1, m
c   1000 write(iunit) (a(i,j), j=1, n)
c  in an ftn77  compiled program
c
      common /inout/ inp,iout,inf,info
      dimension a(mdim,*)
      dimension         lookup(4,4)
      data lookup /1,2,3,4,5,6,7,8,
     $             4,3,2,1,8,7,6,5/
c
c  Transform the number of clockwise rotations to one of  0,1,2,3
c
      klick=mod(mod(klack,4) + 4, 4)
c
c  All possible orientations reduce to one of the 8 below.  lookup
c  has this coded into it
c
      igo=lookup(jsign, klick+1)
      goto (1010, 1020, 1030, 1040, 1050, 1060, 1070, 1080), igo
c
c  No transformations.  Read  a  as written
 1010 do 1015 i=1, m
 1015 read (inf, end=2000) (a(i,j), j=1,n)
      return
c
c  v - reflect in vertical plane
 1020 do 1025 i=m,1,-1
 1025 read (inf, end=2000) (a(i,j), j=1,n)
      return
c
c  h - reflect in horizontal plane
 1030 do 1035 i=1,m
 1035 read (inf, end=2000) (a(i,j), j=n,1,-1)
      return
c
c  vh = r**2  - reflect in v and h =  2 clockwise 90 deg turns
 1040 do 1045 i=m,1,-1
 1045 read (inf, end=2000) (a(i,j), j=n,1,-1)
      return
c
c  r - 1 clockwise 90 deg turn
 1050 do 1055 j=n,1,-1
 1055 read (inf, end=2000) (a(i,j), i=1,m)
      return
c
c  rv - v reflection followed by r
 1060 do 1065 j=n,1,-1
 1065 read (inf, end=2000) (a(i,j), i=m,1,-1)
      return
c
c  rh - h reflection followed by r = transpose of a
 1070 do 1075 j=1,n
 1075 read (inf, end=2000) (a(i,j), i=1,m)
      return
c  r**3=r**-1 - 3 clockwise 90 deg turns = one anticlockwise
 1080 do 1085 j=1,n
 1085 read (inf, end=2000) (a(i,j), i=m,1,-1)
      return
c
c  Error section.  Improper dimensions have been supplied
c
 2000 write(*,'(a)') ' ','>>> Unexpected endfile encountered',
     $'Improper array size has been given -  File rewound'
      m=0
      rewind inf
      return
      end
c______________________________________________________________
      subroutine matraf(mdim, m, n, a, klack, jsign)
c$$$$ calls nothing
c  Routine actually to read formatted diskfile with * format
c  it is assumed the file was written by rows with a statement like
c     write(iunit, format) ((a(i,j), j=1, n), i=1, m)
c
      common /inout/ inp,iout,inf,info
      dimension a(mdim,*)
      dimension         lookup(4,4)
      data lookup /1,2,3,4,5,6,7,8,
     $             4,3,2,1,8,7,6,5/
c
c  Transform the number of clockwise rotations to one of  0,1,2,3
c
      klick=mod(mod(klack,4) + 4, 4)
c
c  All possible orientations reduce to one of the 8 below.  lookup
c  has this coded into it
c
      igo=lookup(jsign, klick+1)
c
      goto (1010, 1020, 1030, 1040, 1050, 1060, 1070, 1080), igo
c
c  No transformations.  Read  a  as written
 1010 read (inf,*, err=2010,end=2000) ((a(i,j), j=1,n), i=1, m)
      return
c
c  v - reflect in vertical plane
 1020 read (inf,*, err=2010,end=2000) ((a(i,j), j=1,n), i=m,1,-1)
      return
c
c  h - reflect in horizontal plane
 1030 read (inf,*, err=2010,end=2000) ((a(i,j), j=n,1,-1), i=1,m)
      return
c
c  vh = r**2  - reflect in v and h =  2 clockwise 90 deg turns
 1040 read (inf,*, err=2010,end=2000) ((a(i,j), j=n,1,-1), i=m,1,-1)
      return
c
c  r - 1 clockwise 90 deg turn
 1050 read (inf,*, err=2010,end=2000) ((a(i,j), i=1,m), j=n,1,-1)
      return
c
c  rv - v reflection followed by r
 1060 read (inf,*, err=2010,end=2000) ((a(i,j), i=m,1,-1), j=n,1,-1)
      return
c
c  rh - h reflection followed by r = transpose of a
 1070 read (inf,*, err=2010,end=2000) ((a(i,j), i=1,m), j=1,n)
      return
c  r**3=r**-1 - 3 clockwise 90 deg turns = one anticlockwise
 1080 read (inf,*, err=2010,end=2000) ((a(i,j), i=m,1,-1), j=1,n)
      return
c
c  Error section.  Improper dimensions have been supplied
c
 2000 write(*,'(a)') ' ','>>> Unexpected endfile encountered',
     $'Improper array size has been given -  File rewound'
      m=0
      if (inf.ne.inp) rewind inf
      return
c
 2010 write(*,'(a)') ' ','>>> Error while reading array data',
     $'File rewound'
      m=0
      if (inf.ne.inp) rewind inf
      return
      end
c______________________________________________________________
      subroutine assess(fmin, fmax, step, nlev, levmx, zlevs)
c$$$$ calls nothing
c  Finds relatively nice subdivisions on interval (fmin, fmax)
c  Returns number of levels in  nlev and array of
c  values in zlevs.
c  If  step > 0 uses this interval to construct the subdivisions
c  but modifies interval to ensure nlev <= levmx
c  step < 0 signals how many contour levels in default scheme.
      dimension zlevs(*),xnice(10)
      data mlines/11/
      data xnice/1,2,2.5,5,5,5,5,10,10,10/
c
      f1=min(fmin, fmax)
      f2=max(fmax, fmin)
      if (f1 .eq. f2) f2=1.0 + 1.5* abs(f1)
c
      numlev=mlines
      if (step.lt. -mlines) numlev=min(-nint(step), levmx)
c
c  Make up a nice uniform contour interval
      if (step.le. 0.0) then
        dzlog = log10((f2-f1)/numlev) + 500.0
        near=nint(10.0**mod(dzlog, 1.0))
        dz = xnice(near)*10.0**(int(dzlog) - 500)
c  Use step provided
      else
        dz = step
      endif
c 
c  Build level lines
      if (f2-f1 .ge. levmx*dz) dz=(f2-f1+dz)/(levmx-1)
      k1 = nint(f1/dz)
      k2 = nint(f2/dz)
      nlev=k2-k1+1
      do 1500 k=k1, k2
        zlevs(k-k1+1)=k*dz
 1500 continue
      return
      end
c______________________________________________________________
      subroutine group(string, k1,k2)
c$$$$ calls nothing
c  Given a string representing repeated application on an array of
c  the actions:
c    R rotate 90 deg clockwise; V reflect in a verical plane
c    H reflect in a horizontal plane; T transpose,
c  finds an integer representing one of the states of the array
c    1=unaltered; 2=V; 3=H; 4=RR; 5=R; 6=T; 7=RH; 8=RRR
c
c  This in turn is decoded into [m n r] form in k1, k2
c
c  The order of operations is as one reads - thus RTV means
c  first rotate array once, then transpose the result, then reflect the
c  result in a vertical plane
      dimension multab(8,8), map(8,2)
      character *(*) string
c  This is the dihedral D2 group multiplication table
      data ((multab(i,j),j=1,8),i=1,8)/
     $1,2,3,4,5,6,7,8,  2,1,4,3,7,8,5,6,  3,4,1,2,6,5,8,7,
     $4,3,2,1,8,7,6,5,  5,6,7,8,4,3,2,1,  6,5,8,7,2,1,4,3,
     $7,8,5,6,3,4,1,2,  8,7,6,5,1,2,3,4/
      data ((map(i,j),j=1,2),i=1,8)/
     $1,1, -1,1, -1,-1, 1,-1, 1,2, -1,2, -1,0, 1,0/
c
      igroup=1
      l=len(string)
      do 1500 i=1, l
        k=index('.VH.RT',string(i:i))
        if(k .ge. 2) igroup=multab(igroup,k)
 1500 continue
c
c  Map orientation into [m, n, r] instructions
      k1=map(igroup,1)
      k2=map(igroup,2)
      return
      end
c________________________________________________________________
*======================= Array input routines END
*====================================
*====================================
*  Map generator and page manager BEGIN ================
*
      subroutine major(mdim, m, n, a)
c$$$$ calls assess bar border clear getchr getone getval letter
c$$$$ trima linax minmax notes paint table sort xyline xysym encaps
c
c  The major routine for producing color map of the array  a.
c  m, n  are the internal array dimensions.
c  mdim  is the row dimension as  a  is stored (=m+2)
c
c  The routine positions the plot on the paper, opens the plotfile,
c  invokes the painting routine and positions the origin, puts
c  the color key at the right and writes the plot label
c
c  Postscript prolog written here also
c
      parameter (levmx=100)
      dimension a(mdim, *), zlevs(levmx), xyorg(2),xynew(2)
      common /xtreme/ bord,amin,amax,alim(3)
      parameter (lpmx=2500, knmx=50 000)
      integer*2 kxy
      common /arcs/ loops,istart(lpmx),area(lpmx),lev(lpmx),kxy(2,knmx)
      common /trace/ level,em,en,xfac,yfac
      common /noir/ scutum,kolors,art(4,levmx),kart(levmx),loutln(levmx)
      common /pallet/ many,npal(10),hue(28,10),sat(28,10),bri(28,10)
      common /inout/ inp,iout,inf,info
      common /local/ color(1),black(2)
      common /outfil/ output
c
      dimension hw(3), xyax(4), zlim(2), rlim(2)
      character*80 output, label, txes, onoff*2,boff*2
      save kase, hw, xyorg, xynew, hlet, legacy, zlevs
      data kase/0/,  hw/2*0, 0.08333/, nlev/0/, legacy/0/
      data xyorg/0,1000/, xynew/1,1.5/
      data onoff/'of'/, boff/'of'/
      data hlet/0.0/
c
c  No array has been provided. Abandon ship
c
      if (m .eq. 0) then
         write(*, '(/a)')'>>> Nothing to plot - read an array'
         return
      endif
c
c  Get a substitute color table for palette 1 if desired
      call table
c
      kase=1 + kase
      loops=0
      istart(1)=1
c  First time through - open plotfile and set landscape/portrait
      if (kase .eq. 1) then
        call getchr('landscape',  label, land)
        output='mypost'
        call getchr('output', output, nf)
        open(unit=9, file=output)
        write(9,'(a)') '%!PS PostScript file',
     $  '% Generated by color',
     $  '%initmatrix','0.072 0.072 scale  % Coords are 1/1000 in'
        if (land .ge. 0) write(9,'(a)')
     $  '6750 0 translate 90 rotate    % anticlockwise',
     $  '%0 10750 translate -90 rotate % clockwise'
        write(9,'(a)') '200 200 translate',
     $  '/m {newpath moveto} def','/r {rlineto}def'
        dxy=200.0/1000
        call encaps(1, dxy,dxy, 0.0,0.0)
      endif
c
c   Determine extremes of function to be plotted; if interval has
c   been set explicitly retain the values
      call minmax(m+2, m ,n, a(2,2), amin, amax)
c
c
c  Set level lines and colors; some contortions needed to inherit
c  colors set by a previous 'above' or 'level' command
c
c  Check whether to use colors previously set 
      if (legacy.gt. 0) then
        call getone('levels', dum,   ilev)
        call getone('interval', dum, iint)
        call getone('above',dum,     nart)
        if (max(iint,ilev,nart).lt. 0) goto 2000
      endif
c
c  Transitions and colors set by "above" commands
      call getval('above', art(1,1), 4, nart)
      call  clear('above')
c
      if (nart .gt. 0) then
        legacy=1
        zlevs(1)=art(1,1)
        do 1100 k=2, levmx
          call getval('above', art(1,k), 4, nart)
          call  clear('above')
          zlevs(k)=art(1,k)
          if (nart .lt. 0) goto 1110
 1100   continue
        write(*,'(a)') '>>> Too many levels requested. ',
     $  'If necessary increase parameter LEVMX and recompile'
        k=levmx + 1
 1110   kolors=max(2, k-1)
        call minmax(kolors, kolors, 1, zlevs, zlim(1), zlim(2))
        kolors=kolors+1
        zlevs(kolors)=zlim(1) - 0.07*(zlim(2)-zlim(1))
c
        call sort(kolors, zlevs, kart)
        nlev=kolors
c  Extrapolate a shade for area below lowest specified level
        k1=kart(1)
        k2=kart(2)
        k3=kart(3)
        art(1,k1)= zlevs(1) 
        art(2,k1)= mod(1.5*art(2,k2) - 0.5*art(2,k3), 1.0)
        art(3,k1)= max(min(1.5*art(3,k2) - 0.5*art(3,k3), 1.0), 0.0)
        art(4,k1)= max(min(1.5*art(4,k2) - 0.5*art(4,k3), 1.0), 0.0)
        goto 2000
      endif
c
c
c  Transitions set by 'interval' command: contours are evenly spaced
      call getval('interval', alim,  3, nlim)
c
c
c  When kolors=0, only transitions have been set, not colors
      kolors=0
      dz=1.0
      if (mod(nlim,2).eq. 1) dz=alim(nlim)
c  Detect various errors in interval specification
      if (nlim.ge.2 .and. alim(2).le.alim(1)) nlim=0
      if (dz.le. 0.0) nlim=0
c  Interval not provided: Default to data extremes; step<0 special case
      if (nlim.le. 0) then
        step=0.0
        call trima(m+2,m,n, a, zlim)
c  Only interval step > 0 provided
      elseif (nlim.eq. 1) then
        legacy=2
        step=alim(1)
        call trima(m+2,m,n, a, rlim)
        zlim(1)=rlim(1)-mod(rlim(1),step) + min(0.0,sign(step,rlim(1)))
        zlim(2)=max(rlim(2), rlim(1)+step)
c  Interval provided, but no step
      elseif (nlim.eq. 2) then
        legacy=2
        step=0.0
        zlim(1)=alim(1)
        zlim(2)=alim(2)
c  Interval and step > 0 provided
      elseif (nlim.eq. 3) then
        legacy=2
        step=alim(3)
        zlim(1)=alim(1) - step
        zlim(2)=alim(2) 
      endif
      if (info .ge. 1) write(*,'(/(a,1p,2g12.4))')
     $' Array values fall within   ',amin,amax,
     $' Contouring performed within',zlim,  ' '
c
c  Assign new levels from 'interval' command
      if (nlim.gt. 0 .or. legacy.eq. 0) then
        call assess(zlim(1), zlim(2), step, nlev, levmx, zlevs)
c  Adjust sub-bottom interval to be same length as that of above-top
        if (nlev.ge. 3) zlevs(1)=zlevs(2) - 0.075*(zlevs(nlev)-zlevs(2))
      endif
c
c
c  Transitions set by 'levels' command 
      call getval('levels', zlevs, levmx, nlevl)
      call  clear('levels')
c
      if (nlevl.gt. 0) then
c  There is a "levels" command
        legacy=3
c  Get any further levels
        do 1500 k=1,4
          call getval('levels', zlevs(nlevl+1), levmx-nlevl, kf)
          call  clear('levels')
          if (kf .le. 0) goto 1510
          nlevl=nlevl + kf
 1500   continue
 1510   nlev=min(nlevl, levmx)
        if (nlev .eq. 1) then
          dist=max(abs(zlevs(1)-amin), abs(zlevs(1)-amax))
          zlevs(2)=zlevs(1) - 1.1*dist
          zlevs(3)=zlevs(1) + 1.1*dist
          nlev=3
        endif
c  Add a sub-bottom level
        call minmax(nlev, nlev, 1, zlevs, zlim(1), zlim(2))
        nlev=nlev+1
        zlevs(nlev)=zlim(1) - 0.07*(zlim(2)-zlim(1))
        zlim(1)=zlevs(nlev)
        call sort(nlev, zlevs, kart)
      endif
c
c
c  Transition levels have been set or inherited from previous plot
c  Choose colors from a standard palette
      call getval('palette', color, 3, nblck)
c  Check alternative spelling "pallet"
      if (nblck.lt. 0) call getval('pallet', color, 3, nblck)
      if (nblck.gt. 0 .or. kolors.eq. 0) then
        kolors=nlev
        match=nint(color(1))
        if (info.ge. 1) write(*,'(a,i3)')
     $    ' Colors drawn from standard palette ',match
c
c  Greytone (match=0): get shade parameters
        if (match.le. 0) then
          bfact=sign(1.0, black(1))/((zlim(2)-zlim(1)) + 1.0e-8)
          offst=0.5 - sign(0.5, black(1))
          ablack=abs(black(1))
          do 1600 k=1, nlev
            art(1,k)=zlevs(k)
            art(3,k)=0.0
            art(4,k)=1.0-ablack* abs(offst +
     $           bfact*(zlevs(k)-zlim(1)))**abs(black(2))
 1600     continue
        else
c 
c  Canned palette HSB values obtained from common /artist/
c  If number of colors required (nlev) < number available (ncav),
c  pick a color from palette, otherwise interpolate
          mtch=min(many, match)
          ncav=npal(mtch)
          beta=(ncav - 1.0001)/(nlev-1)
          do 1700 k=1, nlev
            art(1,k)=zlevs(k)
            alf=1.0 + beta*(k-1)
            i=alf
            if (nlev .ge. ncav) then
              gam=alf-i
              art(2,k)=(1.0-gam)*hue(i,mtch) + gam*hue(i+1,mtch)
              art(3,k)=(1.0-gam)*sat(i,mtch) + gam*sat(i+1,mtch)
              art(4,k)=(1.0-gam)*bri(i,mtch) + gam*bri(i+1,mtch)
            else
              kin=1 + nint((k-1)*(ncav-2.0)/(nlev-1.0))
              if (nlev .eq. ncav) kin=k
              art(2,k)=hue(kin+1,mtch)
              art(3,k)=sat(kin+1,mtch)
              art(4,k)=bri(kin+1,mtch)
            endif
 1700     continue
        endif
        call sort(nlev, zlevs, kart)
      endif
      kolors=nlev
c
c
c
c  Transitions and levels have been set
 2000 if (info .ge. 1) then
        write(*,'(a)')' ',
     $'  These levels and color codes have been set',
     $'     level          hue   staturation brightness'
        do 2050 kl=1, kolors
          k=kart(kl)
          write(*,'(i2,g12.4,3f10.3)') kl, (art(j,k),j=1,4)
 2050   continue
      endif
c
c
c  Fill artificial border with a small value
      bord=min(2.0*amin - amax, 2.0*zlevs(1) - zlevs(nlev))
      do 2100 i=1, m+2
        a(i,1)=  bord
        a(i,n+2)=bord
 2100 continue
      do 2200 j=1, n+2
        a(1,j)=  bord
        a(m+2,j)=bord
 2200 continue
c
c  Decide on the size of the map; it is proportional to the dimensions
c  unless otherwise requested
      call getone('height', hw(2), nf)
      call getone('width',  hw(1), nf)
      w=hw(1)
      h=hw(2)
      if (h .le. 0 .and. w .eq. 0) then
         if (m .ge. n) w=6.0
         if (m .lt. n) h=6.0
      endif
      if (w .le. 0) w=h*(m-1)/(n-1.0)
      if (h .le. 0) h=w*(n-1)/(m-1.0)
c
c  Set the origin
c
c  If the origin of this plot is explicitly set, move it.
c  External units are 1/1000 inch for compatibility with pscrip.
c  Program units are inches
      c=1000
c
      call getval('plot', xyorg, 2, nf)
      if (nf .eq. 2) then
        write(9,'(2f9.2,a)') c*xyorg(1),c*xyorg(2),' translate'
        call encaps(1, xyorg(1),xyorg(2), 0.0,0.0)
      else
c  No explicit origin requested.  Set standard position
        write(9,'(2f9.2,a)') c*xynew(1),c*xynew(2),' translate'
        call encaps(1, xynew(1),xynew(2), 0.0,0.0)
      endif
      xynew(1)=0.0
      xynew(2)=h+1.0
      call encaps(2, 0.0,0.0, w,h)
c
c  Paint the array
      if (info .ge. 1) write(*,'(/a)') ' Painting underway ...'
c
      call paint(m, n, a, zlevs, w, h)
c
      if (info .ge. 1) write(*,'(/a,i5)')
     $' Number of closed contours in map ',loops
c
c  Add extra xy-lines
      call xyline(w, h)
c
c  Add symbols
      call xysym(w, h)
c
c  Draw a border around the map or mask out around inscribed ellipse
      call getchr('border', boff, nbor)
      if (nbor .eq. 0) boff='on'
      if (boff .eq. 'ma') then
        call border(10.0, 3)
        call border(10.0, 2)
      elseif (boff .eq. 'on') then
        map=1
        call getone('mapping', amap, nap)
        if (nap .eq. 1 .and. amap.eq.2) map=2
        call border(10.0, map)
      endif
c
c  Add annotated axes
      naxes=0
      if (hlet .le. 0.0) hlet=int(0.5*h + 9.0)/72.0
      call getone('size', hlet, ignor)
      call getval('axes', xyax, 4, nax)
      if (nax .eq. 4) then
c  Cancel axes if they are hidden or of zero length
        call getchr('axes', txes, nax)
        if (xyax(1).eq.xyax(2) .or. xyax(3).eq.xyax(4) .or.
     $     index(txes, 'hid').gt.0) goto 3000
        naxes=2
        call linax(hlet, (em-1)*xfac/1000, xyax(1), 1)
        call linax(hlet, (en-1)*yfac/1000, xyax(3), 2)
      endif
 3000 continue
c
c  Plot label if there is one
      call getchr('label', label, nf)
c
      if (nf.ge. 1) then
        if (info .ge. 1) write(*,'(/2a)') ' Label: ',label(1:72)
        call letter(0.0, -(naxes+2)*hlet, hlet, 0.0, label)
      endif
c
c
c  Draw a color bar to right of main graph if not suppressed
      call getchr('nobar', onoff, nbar)
      if (nbar.eq. 0) onoff='on'
      if (onoff .eq. 'of') then
        write(9,'(2f9.1,a)') c*(w+hlet), 0.0, ' translate'
        call encaps(1, w+hlet, 0.0, 0.0,0.0)
        call bar(min(5.0, h), hlet)
        write(9,'(2f9.1,a)') -c*(w+hlet), 0.0, ' translate'
        call encaps(1, -w-hlet, 0.0, 0.0,0.0)
      endif
c
c  Plot the notes last
      call notes(hlet)
      return
      end
c______________________________________________________________
      subroutine paint(m, n, a, zlevs, w, h)
c$$$$ calls dumpit getval getchr glevel sort tensor
c  Jacket routine for performing basic contouring function.
c  As always,  a  is the array with original dimensions m, n.
c  zlevs is an array of contour levels; the list must be
c    in ascending numerical order.  Routine quits by detecting
c    a decrease.
c  h, w  are height and width in inches of finished map
      save outln
      parameter (lpmx=2500, knmx=50 000)
      dimension a(*),zlevs(*),basx(103),basy(103),key(lpmx)
      dimension outln(lpmx), heavy(3)
      character*80 txt
      integer*2 kxy
      common /perim/ ltrue(lpmx)
      common /arcs/ loops,istart(lpmx),area(lpmx),lev(lpmx),kxy(2,knmx)
      common /trace/ level,em,en,xfac,yfac
      common /xtreme/ bord,amin,amax,alim(3)
      parameter (levmx=100)
      common /noir/ scutum,kolors,art(4,levmx),kart(levmx),loutln(levmx)
      data txt/' '/
      data heavy/4, 4, 12/
c
      diff(x,y)=max(abs(x-y)-(abs(x)+abs(y))*1e-5, 0.0)
c
c  Inquire about outlines
      noutln=-99
      call getchr('outline', txt, ignor)
      if (index(txt,'off') .eq. 0)
     $  call getval('outline', outln, levmx, noutln)
       do 1050 k=1, kolors
         loutln(k)=abs(noutln)
 1050  continue
c
c  loutln=0 => outline that level
      do 1150 k=1,  kolors
       do 1100 i=1, noutln
         if (diff(zlevs(k), outln(i)) .eq.0.0) loutln(k)=0
 1100   continue
 1150 continue
c
c  Inquire about dashed outlines - marker in loutln set to -1
      call getval('dashed', outln, levmx, noutln)
c  loutln= -1 => dash outline that level
      do 1250 k=1,  kolors
       if (noutln .eq. 0) loutln(k)=-1
       do 1200 i=1, noutln
         if (diff(zlevs(k), outln(i)) .eq.0.0) loutln(k)=-1
 1200   continue
 1250 continue
c  Inquire about heavy outlines - marker in loutln set to -2
      call getval('heavy', outln, levmx, noutln)
      do 1350 k=1,  kolors
       if (noutln .eq. 0) loutln(k)=-2
       do 1300 i=1, noutln
         if (diff(zlevs(k), outln(i)) .eq.0.0) loutln(k)=-2
 1300   continue
 1350 continue
c  Get line weights for outlines
      call getval('weights', heavy, 3, ignor)
c
c  Set scale factors which are single grid steps in x and y
c  and dimensions for dumpit
      c=1000
      xfac=c*w/(m-1)
      yfac=c*h/(n-1)
      em=m
      en=n
c
c  Establish map transformations if any
      call tensor
c
c  Offset the map by one cell width for border
      write(9,'(2f10.1,a)') -xfac,-yfac,' translate'
c
c  Load basx, basy with a bounding rectangle;
c  y values are closely spaced for Aitoff projection
      do 1400 i=1, 51
        basx(i)=em + 1.0
        basy(i)=(en + 1.0)*(1.0 - 0.02*(i-1))
        basx(i+51)=0.0
        basy(103-i)=basy(i)
 1400 continue
c  Determine the base color
      do 1500 level=0, kolors-1
        if (zlevs(level+1) .gt. bord) goto 1510
 1500 continue
c  Fill rectangle with  base color
 1510 level=1
      call dumpit(zlevs(1), 102, basx,basy)
      area(1)=1.5*area(1)
c
c  Set 0.10 inches to be the maximum step on any contour:
      ds=0.10*c/min(xfac, yfac)
c  Run through contour levels, saving paths in  kxy.
c  Save time and memory by skipping useless levels
      do 1600 level=2, kolors
        zlv=zlevs(level)
        if (zlv .ge. amax) goto 1610
        if (zlv .gt. bord) then
          call glevel(m+2, m+2, n+2, a, ds, zlv)
        endif
 1600 continue
 1610 continue
c
c  Sort the contours path into decreasing order by area.
c  Since sort makes ascending lists, numerically area < 0
      call sort(loops, area, key)
c
c  Output the paths in order of decreasing area
      do 2800 lk=1, loops
        l=key(lk)
        level=lev(l)
        nxy=istart(l+1)-istart(l)
        j1=istart(l)   + 1
        j2=istart(l+1) - 1
        ka=kart(level)
        lout=loutln(ltrue(l))
c  Postscript code written to unit 9; each loop annotated with a header
        write(9,'(a,i5,a,i3,a,f9.1,a,g12.4/3f9.3,a)') '%% nxy ',nxy,
     $  '  color number',level,'  area', -area(lk),
     $  '  true level',zlevs(ltrue(l)),
     $  art(2,ka),art(3,ka),art(4,ka),' sethsbcolor'
        write(9,'(2i7,a)') kxy(1,j1-1),kxy(2,j1-1),' m'
        write(9,*) (kxy(1,j) - kxy(1,j-1),
     $              kxy(2,j) - kxy(2,j-1), ' r',j=j1,j2)
        if (lout .gt. 0) then
          write(9,'(a)') 'closepath fill'
c
c  Outline region with a black line
        elseif (lout.eq. 0 .or. lout.eq. -2) then
          write(9,'(a,f8.0)')'closepath gsave fill grestore 0 setgray',
     $    heavy(1-lout),' setlinewidth stroke 4 setlinewidth'
        elseif (lout .eq. -1) then
          write(9,'(2a,f8.0,a)') 'closepath gsave fill ',
     $    'grestore 0 setgray', heavy(2), ' setlinewidth ',
     $    '[50 50] 0 setdash stroke [] 0 setdash 4 setlinewidth'
c
*  A more subtle outliner - hardcopy not up to it though!
*         kb=ltrue(l)
*         write(9,'(3f7.3,a)') art(2,kb+1),art(3,kb+1),
*    $    art(4,kb+1),' sethsbcolor stroke'
*  An alternate that is pretty for grays
*    $    max(0.0,art(4,kb)-0.3),' setgray stroke'
        endif
 2800 continue
c
c  Restore origin and color to black
      write(9,'(2f10.1,a)')  xfac,yfac,' translate 0 setgray'
      return
      end
c____________________________________________________________________
      subroutine border(weight, kind)
c$$$$ calls nothing
c  Draws a border around the finished map
c  weight  is line width in 1/1000th inches
      common /trace/ level,em,en,xfac,yfac
      common /atlas/ map,mtab,ntab,xtab(201),ytab(201)
      data f/1.004/
c
c  Define sizes
      a=0.5*(em - 1.0)*xfac
      b=0.5*(en - 1.0)*yfac
      kolor=(kind-1)/2
      wt=abs(weight)
c
c  Choose appropriate shape depending on the map
      goto (1000, 2000, 1000), kind
c
c  Map has rectangular border (white for masking)
 1000 write(9,*) '% Rectangular border'
      write(9,*) 'gsave ',kolor,' setgray ',wt,' setlinewidth'
      write(9,*) 0,0,' m',2*a,0,' r',0,2*b,' r',-2*a,0,' r'
      write(9,*) 'closepath stroke grestore'
      if (kind .eq. 3) goto 3000
      return
c
c  Elliptical border - drawn by the PostScript arc command
 2000 write(9,*) '% Elliptical border'
      write(9,*) 'gsave',a, b, ' translate 1 ',b/a,' scale '
      write(9,*) 'newpath 0 0',a,' 0 360 arc'
      write(9,*) wt,' setlinewidth 0 setgray stroke grestore'
      return
c
c  Mask out stuff between the inscribed ellipse and the rectangle
 3000 write (9,*) 'gsave', a,b,' translate',a*f,b*f,' scale'
      write (9,*) 'newpath 0 0 ',1.0/f,
     $            ' 0 360 arc 0 -1 r -2 0 r 0 2 r 2 0 r'
      write (9,*) '0 -1 r closepath 1 setgray fill grestore'
      return
c
      end
c______________________________________________________________________
      subroutine bar(height, hlett)
c$$$$ calls getval getchr letter mabel
c  Draws the color key to the right of the map.  This involves
c  numerical annotation at the proper place done by mabel
c
      character*80 label
      common /xtreme/ bord,amin,amax,alim(3)
      parameter (levmx=100)
      common /noir/ scutum,kolors,art(4,levmx),kart(levmx),loutln(levmx)
c
      common /inout/ inp,iout,inf,info
      dimension heavy(3), rung(levmx+1)
      data heavy/4, 4, 12/, label/' '/
c
c
      call getval('weights', heavy, 3, ignor)
      call getchr('size', label, ignor)
      strict=min(index(label, 's'), 1)
      c=1000
      hlet=hlett
      art1=art(1,kart(1))
      kk=kart(kolors)
      cons=0.94*c*height/(1.0e-18 + (art(1,kk) - art1))
c  Paint the highest color; just paint the rectangle
c  in the proper color and fill
      dx=2.0*c*hlet
      dz=c*height
      write(9,*) '% Color bar'
      write(9,'(4(a,f8.0))') '0 0 m',dx,' 0 r 0',dz,' r',-dx,' 0 r'
      write(9,'(a,3f10.3,a)') 'closepath ',
     $art(2,kk),art(3,kk),art(4,kk),' sethsbcolor fill'
c
c  Determine minimum contour spacing, and set initial spacing of
c  lines in the bar
      drmin=1.0e6
      do 1100 k=1, kolors
        rung(k)=cons*(art(1,kart(k)) - art1)
        if (k.ge. 2) drmin=min((rung(k)-rung(k-1))/c, drmin)
 1100 continue
      rung(kolors+1)=dz
      deven=(rung(kolors)-rung(2))/(kolors-2)
c
c  Check for crowding; then adjust height of numbers and, if necessary, 
c  stretch scale
      shrnk=min(0.83*deven/(c*hlett), 1.0)
      if (drmin .lt. 1.2*hlet) then
        hlet=shrnk*hlett
        if (drmin .ge. 1.2*hlet) goto 1200
        gamma=deven
        alf=(1.2*hlet*c-gamma)/(drmin*c-gamma)
c  Distort intervals on bar to fit numbers in without overlap
        do 1150 k=2, kolors
          rung(k)=alf*rung(k) + (1.0-alf)*(gamma*(k-2)+rung(2))
 1150   continue
      endif
 1200 continue
c
c  Fill the colors and the outlines
      do 1400 kl=kolors, 2, -1
        ka=kart(kl)
        dz=rung(kl)
        write(9,'(4(a,f8.0))') '0 0 m',dx,' 0 r 0',dz,' r',-dx,' 0 r'
        ka=kart(kl-1)
        write(9,'(a,3f10.3,a)') 'closepath ',
     $  art(2,ka),art(3,ka),art(4,ka),' sethsbcolor fill'
c  Take care of outline marker in color bar
        klout=loutln(kl)
        if (klout.eq. 0 .or. klout .eq. -2) write(9,'(a,f8.0)') 
     $     '0 ',dz,' m',dx,' 0 r 0 setgray ',heavy(1-klout),
     $     ' setlinewidth stroke 4 setlinewidth'
        if (loutln(kl) .eq.-1) write(9,'(a,f8.0)')
     $     '0 ',dz,' m',dx,' 0 r ', heavy(2),
     $     ' setlinewidth 0 setgray [50 50] 0 setdash stroke' //
     $     ' [] 0 setdash  4 setlinewidth'
 1400 continue
c
c  Surround the bar with a border
      write(9,'(4(a,f8.0))')
     $'0 0 m',dx,' 0 r 0',height*c,' r',-dx,' 0 r'
      write(9,'(a)') '0 setgray closepath stroke'
      write(9,'(a)') '0 setgray'
c
c  Restore original letter height if sizes are strict 
      if (strict.eq. 1.0) hlet=hlett
      iskip=1.0 + strict*0.9/shrnk
c
c  If numbers are too large or small inmagnitude, scale them
      larg=log10(max(abs(art(1,kk)), abs(art1)))
      if (larg .lt. 0) larg=larg - 1
      scale=1.0
      dx=dx + c*hlet
c  Write the scale factor in exponential notation above the list
      if (larg .ge. 5 .or. larg .le. -2) then
        scale=10.0**larg
        dz=cons*(art(1,kk) - art1) + 0.75*c*hlet
        call letter(dx/c, dz/c, hlet, 0.0, 'x10')
        call mabel(real(larg), klen, label)
        call letter(dx/c+1.85*hlet, dz/c+0.5*hlet, 0.75*hlet, 0.0,
     $              label(1:klen))
      endif
c
c  Write numerical annotation to the right of bar
      do 1500 k=2, kolors, iskip
        ka=kart(k)
        call mabel(art(1,ka)/scale, klen, label)
        call letter(dx/c, rung(k)/c-0.5*hlet, hlet, 0.0, label(1:klen))
 1500 continue
      if (iskip.gt. 1 .and. drmin*c .lt. 0.9*deven) write(*,'(a)')
     $'>>> Color bar information loss: unevenly spaced unlabeled levels'
c
      return
      end
c______________________________________________________________________
      subroutine letter(x, y, h, angle, text)
c$$$$ calls greek encaps
      character*(*) text, bkslsh*1
c  A PostScript version of the famous routine, without 
c  subscript ability.
c  External units are reset to points (1/72 inches);
c  Program units as always are inches
c
      bkslsh=char(92)
      long=len(text)
      do 1100 l=long, 1, -1
        if (text(l:l) .ne. ' ') goto 1101
 1100 continue
      return
 1101 long=l
c
      ipoint=nint(h*72)
      write(9,'(a)')'gsave 15.0 13.89 scale 0 setgray'
c
      write(9,'(a)') '/Helvetica findfont'
      write(9,'(i4,a/2f6.1,a)') ipoint,' scalefont setfont',
     $          67*x, 72*y, ' moveto currentpoint translate'
      if (angle .ne. 0) write(9,'(f8.2,a)') angle,' rotate'
c
      i1=1
 1000 continue
      i2=i1 -1 + index( text(i1:long),  bkslsh)
      if (i2.eq. i1-1) i2=long+1
      if (i2 .gt. i1) write(9,'(3a)') '(',text(i1:i2-1),') show'
      if (i2 .ge. long) goto 1500
      i3=i2 + index( text(i2+1:long), bkslsh)
      if (i3.eq. i2) i3=long+1
      if (i3-i2.gt. 1) call greek(ipoint, text(i2+1:i3-1))
      if (i3.ge. long) goto 1500
c
      i1=i3+1
      goto 1000
c
c
 1500 write(9,'(a)') 'grestore'
c
      el=0.6*l*h
      px=el*cos(0.01745*angle)
      py=el*sin(0.01745*angle)
      call encaps(2,x,y, x+px,y+py)
      call encaps(2,x-py*h/el, y+px*h/el, x+px-py*h/el, y+py+px*h/el)
c
      return
      end
c______________________________________________________________________
      subroutine greek(ipoint, text)
c$$$$ calls nothing
c  Draws a single Greek character; then returns to Roman font.
c  Also interprets \up\ and \down\ as half height vertical motions
      character*(*) text
c
c  Check superscripting
      if (text .eq. 'up') then
        write(9,'(2i5,a)') 0,ipoint/2,' r'
        return
c  Check subscripting
      elseif (text .eq. 'down') then
        write(9,'(2i5,a)') 0,-ipoint/2,' r'
        return
      endif
c
c  Resolve ambiguities in representation of letters.
c  Symbol font maps Roman to Greek with a few oddities handled below
      if (text(1:2).eq.'et') text(1:1)='h'
      if (text(1:2).eq.'ph') text(1:1)='f'
      if (text(1:2).eq.'ps') text(1:1)='y'
      if (text(1:3).eq.'ome') text(1:1)='w'
      if (text(1:2).eq.'Et') text(1:1)='H'
      if (text(1:2).eq.'Ph') text(1:1)='F'
      if (text(1:2).eq.'Ps') text(1:1)='Y'
      if (text(1:3).eq.'Ome') text(1:1)='W'
c
      write (9,*)'/Symbol findfont ',ipoint,' scalefont setfont'
      write (9,'(3a)') '(',text(1:1),') show'
      write (9,*)'/Helvetica findfont ',ipoint,' scalefont setfont'
      return
      end
c_______________________________________________________________________
      subroutine notes(hlet)
c$$$$ calls getchr letter clear
c  Writes a series of notes at the given coordinates (inches) with
c  option letter height, and clockwise angle parameters
      common /inout/ inp,iout,inf,info
      character*75 line
c
      do 1800 k=1, 1000
         call getchr('note', line, nnote)
         call  clear('note')
         if (nnote.lt. 0) then
           if (k.gt.1) write(*,'(/a,i4)')' Notes added to map:',k-1
           return
         endif
         i1=index(line, '(')
         i2=index(line,')')
         if (i1*i2 .eq. 0 .or. i2.lt.i1) goto 1700
         read(line(i1+1:i2-1), *, end=1050, err=1700) x,y,h,angle
         if (h.eq. 0.0) h=hlet
         call letter(x,y, h,angle, line(i2+1:75))
         goto 1800
 1050    read(line(i1+1:i2-1), *, end=1100, err=1700) x,y,h
         if (h.eq. 0.0) h=hlet
         if (abs(h).le. 2) call letter(x,y, h,0.0, line(i2+1:75))
         goto 1800
 1100    read(line(i1+1:i2-1), *, err=1700) x,y
         call letter(x,y, hlet,0.0, line(i2+1:75))
         goto 1800
 1700    write(*,'(a,a)')'>>> Improper note coordinates: ',line
 1800 continue
      return
      end
c_______________________________________________________________________
      subroutine table
c$$$$ calls getchr
c  Reads in a palette from a file, substituting it for
c  the standard palette number 1.  Always reads to EOF
c
      parameter (nk=28)
      character*64 name
      common /inout/ inp,iout,inf,info
      common /pallet/ many,npal(10),hue(nk,10),sat(nk,10),bri(nk,10)
c
c  Look for filename of a color table
      call getchr('table', name, ntb)
      if (ntb .le. 0) return
      open (unit=10, file=name, status='OLD', err=2000)
c
c  Read into the first palette
      do 1100 it=1, nk
        read (10,*,end=1110,err=2100) hue(it,1),sat(it,1),bri(it,1)
 1100 continue
      it=nk + 1
 1110 it=it - 1
      close(unit=10)
      npal(1)=it
      write(*, '(/a,i3,a)')
     $' Palette 1 has been replaced by a new selection of ',it,' colors'
      return
c
c  Errors reported
 2000 write(*, '(a,a/a)')'>>> File ',name,
     $'is nonexistent or you lack permission to read it'
      return
 2100 write(*, '(2a)')'>>> Unreadable numbers in file ', name,
     $'Palette 1 will be a mixture of old and new colors'
      return
      end
c_______________________________________________________________________
      subroutine minmax(mdim, n, k, x, xmin, xmax)
c$$$$ calls nothing
c  Finds minimum and maximum of an array  x  with actual dimensions
c  n, k  held in columns that are mdim long
      dimension x(mdim,*)
      xmin=x(1,1)
      xmax=xmin
      do 1100 i=1, n
        do 1000 j=1, k
          xmin=min(xmin, x(i,j))
          xmax=max(xmax, x(i,j))
 1000   continue
 1100 continue
      return
      end
c______________________________________________________________________
      subroutine trima(mdim, m,n, a, rlim)
c$$$$ calls nothing
c  Finds trimmed extremes, rlim, of a, based on the 1-norm, after
c  trimming extremes over 95% of normal.
      common /xtreme/ bord,amin,amax,alim(3)
      dimension a(mdim,*),rlim(2)
c
c  Scan the array 4 times
      do 1400 ifour=1, 4
        sum=0.0
        newn=0
        do 1100 i=1, m
          do 1000 j=1, n
            goto (100,200,300,300), ifour
 100        sum=sum+a(i,j)
            goto 1000
 200        sum=sum+abs(a(i,j)-abar)
            goto 1000
 300        if (a(i,j).lt. top .and. a(i,j).gt. bot) then
              newn=newn+1
              if (ifour.eq. 3) sum=sum+a(i,j)
              if (ifour.eq. 4) sum=sum+abs(a(i,j)-abar)
            endif
 1000     continue
 1100   continue
        if (ifour.eq. 1)abar=sum/(m*n)
        if (ifour.eq. 3)abar=sum/newn
        if (ifour.eq. 2) then
          bot=abar-2.9*sum/(m*n)
          top=2.0*abar-bot
        endif
 1400 continue
      awidth=sum/newn
      rlim(1)=max(amin, abar-2.8*awidth)
      rlim(2)=min(amax, abar+2.8*awidth)
      return
      end
c_______________________________________________________________________
      subroutine sort(no, x, key)
c$$$$ calls nothing
c  In-place rapid sorter.  Sorts array in place and saves an index to
c  the sort in key
      dimension x(*), key(*)
c
      do 1 i=1,no
 1    key(i)=i
      mo=no
 2    if (mo-15) 21,21,23
 21   if (mo.le.1) return
      mo=2*(mo/4) + 1
      goto 24
 23   mo=2*(mo/8) + 1
 24   ko=no - mo
      jo=1
 25   i=jo
 26   if (x(i) - x(i+mo)) 28,28,27
 27   temp=x(i)
      x(i)=x(i+mo)
      x(i+mo)=temp
      kemp=key(i)
      key(i)=key(i+mo)
      key(i+mo)=kemp
      i=i - mo
      if (i-1) 28,26,26
 28   jo=1 + jo
      if (jo-ko) 25,25,2
      end
c______________________________________________________________________
*==================  Map generator and page manager END
*====================================
      subroutine linax(ht, xlen, xx, ixy)
c$$$$ calls justy letter nicer plot encaps
c  Plots either an x- or a y-axis of length  xlen  inches with tick
c  marks at reasonable places between the limits  xx(1), xx(2).
c  If  xx(1) .gt. xx(2) the numerals decrease along the axis.
c  lettering height is  ht.
c  ixy =1 means  x axis,  ixy=2  means y axis.
c  The axis is drawn starting from the plot origin (x=0, y=0), which is
c  normally at the lower left corner.  This position may be moved by
c  the user in the calling routine
c
      character*40 label,ifmt*20
      dimension xx(2),xa(2),s(4)
c
c
      data ticfac/1.6/, iup,idn /3,2/
c
c  In-line function
      iten(x)=int(500.001 + log10(x)) - 500
c
ce
      if (xlen .eq. 0.0) return
      revers=sign(1.0, xx(2) - xx(1))
      xx(1)=revers*xx(1)
      xx(2)=revers*xx(2)
      htix=0.4*ht
      xa(1)=xx(1)
      xa(2)=ticfac*min(1.0,10.0/xlen)*(xx(2)-xx(1)) + xx(1)
      call nicer(xa, xa, xt)
      if (xx(2)-xx(1).eq.180 .or. xx(2)-xx(1).eq.360) xt=30
c  Set up format for numerical annotation on the axis
      nsiz=iten(max(abs(xx(2)), abs(xx(1))))
      ntix=iten(xt)
      nfld=3 + max(nsiz, -ntix, nsiz-ntix)
      if (nfld.gt.7) goto 1000
c  Width of field less than 8 characters - use an f-format
      ndp=max(0, -ntix)
      write(ifmt, 100) '(f', nfld, '.', ndp, ')'
 100  format(a, i1, a, i1, a)
      if (ndp .eq. 0) nfld=nfld - 1
      goto 1100
c  Width of field more than 7 characters - use a g format
 1000 ndp=nsiz - ntix + 1
      nfld=7 + ndp
      write(ifmt, 110) nfld, ndp
 110  format('(1pg', i2, '.', i2, ')' )
c
c  Reset origin
 1100 call plot(0.0, 0.0, -1)
      xs= xlen/(xx(2) - xx(1))
      n1=xa(1)/xt
      n2=xx(2)/xt + 1.0
      goto (2100, 3100), ixy
c
c  Draw the x axis
 2100 kskip=1.3 + ht*nfld/(xs*xt)
      do 2500 n=n1, n2
        x=n*xt
        xinch=xs*(x - xx(1))
c  Plot next section of axis and the tick on the right
        if (xinch.lt.-0.01 .or. xinch.gt.xlen+0.01) goto 2500
        call plot(xinch, 0.0,  idn)
        call plot(xinch, htix, idn)
c  Write numerical annotation
        if (mod(n, kskip) .ne. 0) goto 2400
        write(label, ifmt) x*revers
        call justy(s, ht, label(1:nfld))
        call letter(xinch-s(2), -1.7*ht, ht, 0.0, label(1:nfld))
c  Move back onto axis with pen up
 2400   call plot(xinch, 0.0, iup)
 2500 continue
      call encaps(2, 0.0, -1.7*ht, 0.0, 0.0)
c  Plot the last little piece of the axis
      call plot(xlen,  0.0, idn)
      call plot(0.0,   0.0, iup)
      xx(1)=revers*xx(1)
      xx(2)=revers*xx(2)
      return
c
c  Draw the y axis
 3100 kskip=1.5 + ht/(xs*xt)
      do 3500 n=n1, n2
        y=n*xt
        yinch=xs*(y - xx(1))
c  Plot next section of axis and the tick on the right
        if (yinch.lt.-0.01 .or. yinch.gt.xlen+0.01) goto 3500
        call plot(0.0,  yinch, idn)
        call plot(htix, yinch, idn)
c  Write numerical annotation
        if (mod(n, kskip) .ne. 0) goto 3400
        write(label, ifmt) y*revers
        call justy(s, ht, label(1:nfld))
        call letter(-ht/2-s(3), yinch-ht/3, ht, 0.0, label(1:nfld))
 3400   call plot(0.0, yinch, iup)
 3500 continue
      call encaps(2, -ht/2-s(3), 0.0, 0.0,0.0)
c  Plot last little piece of axis
      call plot(0.0,  xlen, idn)
      call plot(0.0,  0.0,  iup)
      xx(1)=revers*xx(1)
      xx(2)=revers*xx(2)
      return
      end
c______________________________________________________________
      subroutine nicer(xin, xout, xtick)
c$$$$ calls nothing
c  Routine for scaling intervals and providing tick marks for axes.
c  between 7 and 15 ticks are made, which is suitable for 10in axes.
c    Input parameters
c  xin(1),xin(2)  extremes of variable  x  in its own units
c    output parameters
c  xout(1),xout(2)  adjusted extremes, made to be round numbers.  Note
c    the new interval always covers old one.
c  xtick  distance between tick marks in  x  units (not inches).  This
c    number is always a round number.
      dimension  divs(4),xin(2),xout(2)
c
      data e/1.0e-7/,divs/.1,.2,.5,1.0/
c
ce
      xout(1)=xin(1)
      xout(2)=xin(2)
      if (xout(2).eq.xout(1)) xout(2)=1.0 + 1.1*xout(2)
      plus=1000.0 + log10(xout(2)-xout(1))
      index=1.4969 + 2.857*mod(plus, 1.0)
      units=divs(index)*10.0**(int(plus)-1000)
      bias=(1.+e)*units*aint(1.+max(abs(xout(1)), abs(xout(2)))/units)
      xout(1)=xout(1) - mod(xout(1)+bias, units)
      xout(2)=xout(2) - mod(xout(2)-bias, units)
      if (abs(xout(1)/units) .le. .01) xout(1)=0.0
      if (abs(xout(2)/units) .le. .01) xout(2)=0.0
      xtick=units
      return
      end
c______________________________________________________________________
      subroutine plot(x, y, ipen)
c$$$$ calls nothing
c  Thin version of Calcomp to PostScript pen motion controller
      save moves
      data moves/0/, init/0/
c
      if (init .eq. 0) then
        init=1
*  Terminal moveto removed from definition of c
        write(9,'(a)')'/l {lineto} def',
     $   '/c {currentpoint stroke moveto newpath moveto} def',
     $   '/n {newpath moveto} def'
      endif
      ix=1000.0*x
      iy=1000.0*y
c
c  Extend current polygon
      if (ipen .eq. 2) then
        write(9,*) ix,iy,' l'
        moves=moves + 1
c  Terminate previous line, go to new position ready for new one
      elseif (ipen .eq. 3) then
        if (moves .gt. 0) then
          write(9,*) ix,iy,' c'
          moves=0
        else
          write(9,*) ix,iy,' n'
        endif
c  Reset plot origin and go there
      elseif (ipen .lt. 0) then
         write(9,*) ix,iy,' translate 0 0 moveto'
         moves=0
      elseif (ipen .eq. 6) then
         write(9,*) ix,iy,' lineto closepath fill'
         moves=0
      endif
      return
      end
c_______________________________________________________________________
      subroutine justy(s, ht, text)
c$$$$ calls nothing
c  Thin version of text justification routine
      dimension s(4)
      character *(*) text
c
      w=0.40*ht
      n=len(text)
      do 1100 j=1, n
        if (text(j:j) .ne. ' ') goto 1150
 1100 continue
      s(1)=0.0
      s(2)=0
      s(3)=0
      s(4)=0
      return
 1150 s(1)=(j-1)*w
c
      j1=j+1
      do 1200 j=n, j1, -1
        if (text(j:j) .ne. ' ') goto 1250
 1200 continue
      j=j1
 1250 s(3)=j*w
      s(2)=(s(1)+s(3))/2
      s(4)=n*j
c
      return
      end
c_______________________________________________________________________
      subroutine symbl(x, y, h, n)
c$$$$ calls plot
c  Simple symbol generator matching Calcomp series approximately.
c  Every symbol is either a closed polygon (possibly filled), or
c  a set of evenly spaced rays from the center.
c  0 square; 1 triangle; 2 octagon; 3 diamond; 4 plus; 5 asterisk;
c  6 cross; 7 5-ray; 8 Y base up; 9 pentagon; 10 triangle, base up;
c  11 hexagon; 12 Y; 13 bar; 14 6-ray; 15 dot; 16 heptagon;
c  17 circle; 18 lg circle; 19 sm filled circle; 20 sm filled square;
c  21 sm filled triangle.
c  Code for index n :
c  +-hh aa nn (eg 10 04 04);    + = polygon, - = ray
c  hh/10 = height factor scaled with h for final figure
c  2*pi/aa = initial angle for 1st point
c  nn = number of side of the polygon, or number of rays
c
      dimension kode(22)
      data kode/ 10 08 04, 10 04 03, 08 16 08, 10 01 04,
     $          -10 01 04,-07 04 09,-10 08 04,-09 04 05,
     $          -10 04 03, 10 04 05, 10 12 03, 10 04 06,
     $-08 12 03,-10 04 02,-09 04 06, 02 04 03, 09 04 07,
     $           10 01 16, 20 01 20, 06 01 10, 08 08 04,
     $           10 04 03/
      data twopi/ 6.2832/, thet/0/
c
      m=min(22, max(1, n+1))
      kind=sign(1, kode(m))
      k=abs(kode(m))
      r=h*k*0.5e-5
      phi=twopi/(mod(k,10000)/100)
      nn=mod(k, 100)
c
      call plot(x+r*cos(phi), y+r*sin(phi), 3)
      do  1500 j=1, nn
        if (kind .lt. 0) call plot(x, y, 2)
        thet=phi + j*twopi/real(nn)
        call plot(x+r*cos(thet), y+r*sin(thet), 2)
 1500 continue
      if (n .ge. 19) call plot(x+r*cos(thet), y+r*sin(thet), 6)
      call plot(x, y, 3)
      return
      end
c_______________________________________________________________________
      subroutine xysym(w, h)
c$$$$ calls aitoff axchck getchr getval clear symbl
c  Inputs x-y symbol data from a series of files
c  Establishes scale through the axes command
c
      character*80 list
      parameter (nsymx=5000, ksmx=20)
      save xy,axes
      dimension axes(4),xy(nsymx),hsb(3,8)
      dimension kind(ksmx),symht(ksmx),kolr(ksmx),kbeg(ksmx+1)
      common /inout/ inp,iout,inf,info
      data axes/0,0,0,0/, nax/-99/, kbeg(1)/1/, nsym/0/
      data hsb/0,0,0, 0,0,1, 1,1,1, 0.667,1,1, 0.4,1,0.8,
     $0.15,1,1, 0.1,0.8,1, 0.833,1,1/
c
c
      call axchck(linmap, axes)
c
c  Execute the symbol command up to ksmx times
      do 1200 ksym=1, ksmx
        call getchr('symbol', list, nlist)
        call  clear('symbol')
c  File list exhausted or blank
        if (nlist .lt. 0) then
          if (ksym .gt. 1) goto 1250
          if (nsym .eq. 0) return
          goto  1300
        elseif (nlist .eq. 0) then
          write(*,*)'x-y symbol data erased'
          nsym=0
          return
        endif
c  Decode filename, symbol kind and height and possibly color
        jsplit=index(list, ' ')
        kind(ksym) =9
        symht(ksym)=0.1
        read (list(jsplit:80),*, err=1035) kind(ksym),symht(ksym)
        ksp=index(list,'=') + 1
        kolr(ksym)=1
        if (ksp .gt. 1) kolr(ksym)=1 +
     $     index('blawhiredblugreyelorapur',list(ksp:ksp+2))/3
        goto 1040
 1035   write(*,'(a)')'>>> Cannot read symbol height and/or kind',
     $ '    Arbitrarily assign kind=9, height=0.1'
c
 1040   list(jsplit:80)=' '
        open (unit=11, file=list, status='OLD', err=1045)
        goto 1050
 1045   write(*,'(a,a)')'>>> Nonexistent file: ',list
        kbeg(ksym+1)=kbeg(ksym)
        goto 1200
c
c  Read file to eof
 1050   do 1100 j=kbeg(ksym), nsymx, 2
          read (11,*, end=1150, err=1160) xy(j),xy(j+1)
          if (linmap.gt. 0) call aitoff(xy(j),xy(j+1), xy(j),xy(j+1))
 1100   continue
        j=nsymx
        write(*,'(a,i5,a)') '>>> Too many symbols (>', nsymx/2,')',
     $ '    Excess ignored'
 1150   close(unit=11)
        if (info .ge. 1) write(*,'(2a)')
     $   ' X-y symbol data read from: ',list
        kbeg(ksym+1)=j
        goto 1200
 1160   write(*,'(a,a)')
     $   '>>> Error reading x-y symbol data in file: ',list
        kbeg(ksym+1)=j
        close(unit=11)
c
 1200 continue
      ksym=ksmx+1
      write(*,'(a)')'>>> Too many symbol files',
     $'    Excess ignored'
 1250 nsym=ksym-1
c
c
c  Check if the axes have been defined
 1300 call getval('axes', axes, 4, nax)
c
      if (axes(1).ne.axes(2) .and. axes(3).ne.axes(4)) then
        c=1000
        xfac=c*w/(axes(2) - axes(1))
        yfac=c*h/(axes(4) - axes(3))
        write(9,*)'/l{lineto}def 0 setgray'
      else
        write(*,'(a)')
     $  '>>> The axes command must be provided for symbols',
     $  '    x-y symbol data can not be plotted'
        return
      endif
c
c
c  Write the scaled points to 9
      write (9,*)'% Begin x-y symbol data'
      do 1500 ksym=1, nsym
        write(9,'(3f8.2,2a)') (hsb(j,kolr(ksym)),j=1,3),' sethsbcolor'
        do 1400 j=kbeg(ksym), kbeg(ksym+1)-2, 2
          if (xy(j).ge.axes(1).and.xy(j).le.axes(2)
     $       .and.xy(j+1).ge.axes(3).and.xy(j+1).le. axes(4))
     $    call symbl(xfac*(xy(j)-axes(1))/c, yfac*(xy(j+1)-axes(3))/c,
     $              symht(ksym), kind(ksym))
 1400   continue
 1500 continue
      return
c
      end
c_______________________________________________________________________
      subroutine axchck(linmap, axes)
c$$$$ calls getchr getone getval
c  Returns linmap=1 if data mapping requested, 0 otherwise
c  Also sets axes for mapping 2
      character*80 mapcom
      dimension axes(4)
c
c  Automatically set axes if map 2 (Aitoff) is set
      call getval('axes', axes, 4, nax)
      amap=0
      call getone('mapping', amap, nap)
      if (nax.le. 3 .and. amap.eq. 2.0) then
        axes(1)=-2.0
        axes(2)= 2.0
        axes(3)=-1.0
        axes(4)= 1.0
      endif
c   Check if point and line data are to be mapped too
      call getchr('mapping', mapcom, ignor)
      linmap=0
      if (amap.eq. 2.0) linmap=min(1, index(mapcom,'+'))
      return
      end
c______________________________________________________________________
      subroutine xyline(w, h)
c$$$$ calls aitoff axchck getchr getval clear
c  Inputs then immediately plots xy-line data from a file or
c  series of files; establishes scale through the axes command
c
      character*80 list,mark*2
      parameter (nxymx=17270, mxln=64)
      save xy,axes,nxy,lwt,kolr
      dimension axes(4),xy(nxymx),hsb(3,8)
      dimension kolr(mxln),lin(mxln+1),lwt(mxln)
      common /inout/ inp,iout,inf,info
      data axes/0,0,0,0/, nax/-99/, nxy/0/, lin(1)/1/
      data hsb/0,0,0, 0,0,1, 1,1,1, 0.667,1,1, 0.4,1,0.8,
     $0.15,1,1, 0.1,0.8,1, 0.833,1,1/
c
      call axchck(linmap, axes)
c
c  Execute the line command, repeatedly; read line data
      do 1500 newln=1, mxln
        call getchr('lines', list, nlist)
        call  clear('lines')
c
        if (nlist .lt. 0) then
c  No more lines commands - exit loop
          goto 1520
        elseif (nlist .eq. 0) then
          write(*,*)'x-y line data erased'
          nxy=0
          return
        else
c  Read from the named file;  open a new file
c  Clear previous line data, but continue
          if (newln .eq. 1) then
            if (nxy .gt. 0) write(*,'(a)')
     $     ' Previous x-y line data erased'
            nxy=0
          endif
          jsplit=index(list, ' ')
c  Get line color
          kolr(newln)=1
          ico=index(list(jsplit:80), 'co')
          ieq=index(list(jsplit+ico:80), '=')
          ksp=ico + ieq + jsplit
          if (ico*ieq.gt. 0) kolr(newln)=1 +
     $      index('blawhiredblugreyelorapur',list(ksp:ksp+2))/3
c  Get line weight
          iwt=index(list(jsplit:80), 'wt')
          ieq=index(list(jsplit+iwt:80), '=')
          lwt(newln)=6
          if (iwt*ieq.gt. 0)
     $    read(list(jsplit+iwt+ieq:80),*,end=1040,err=1040) lwt(newln)
 1040     list(jsplit:80)=' '
c  Read from the named file;  open a new file
          open (unit=11, file=list, status='OLD',iostat=is, err=1490)
          do 1050 j=nxy+1, nxymx, 2
            read (11,*, end=1100, err=3100) xy(j),xy(j+1)
            if (linmap.gt. 0) call aitoff(xy(j),xy(j+1), xy(j),xy(j+1))
 1050     continue
          j=nxymx
          write (*,'(a,i5,a)')
     $  '>>> x-y line series truncated to ', nxymx/2,' terms'
 1100     close(unit=11)
          nxy=j-1
          lin(newln+1)=nxy+1
          if (info .ge. 1) write(*,'(/2a/a,i7)')
     $   ' x-y line data read from: ',list(1:index(list, ' ')),
     $   ' Number of points:', (lin(newln+1)-lin(newln))/2
        endif
 1490   if (is.ne. 0)write(*,'(a,a)')'>>> Nonexistent file: ',list
 1500 continue
      write(*,'(a)')'>>> No more x-y line files can be accepted'
 1520 continue
      if (nxy .eq. 0) return
c
c
c  Check if the axes have been defined
      call getval('axes', axes, 4, nax)
c
      if (axes(1).ne.axes(2) .and. axes(3).ne.axes(4)) then
        c=1000
        xfac=c*w/(axes(2) - axes(1))
        yfac=c*h/(axes(4) - axes(3))
        write(9,*)'/l{lineto}def 0 setgray'
      else
        write(*,'(a)')
     $  '>>> The axes command is needed to provide a scale for lines',
     $  '    x-y line data can not be plotted'
        return
      endif
c
c
c  Write scaled lines to unit 9; break line if point goes outside box
      write (9,*)'% Begin x-y line data'
      iseg=1
      do 1800 j=1, nxy, 2
        if (j.eq. 1 .or. j.eq.lin(iseg)) then
          if (j.gt. 1) write(9,*)'stroke'
          write(9,'(3f8.2,a)') (hsb(k,kolr(iseg)),k=1,3),' sethsbcolor'
          write(9,'(i4,a)') max(1, min(150, lwt(iseg))),' setlinewidth'
          iseg=iseg + 1
          mark=' m'
        endif
        if ((xy(j  )-axes(2))*(xy(j  )-axes(1)).gt.0.0 .or.
     $      (xy(j+1)-axes(4))*(xy(j+1)-axes(3)).gt.0.0) then
          if (mark .eq. ' l') write(9,*)'stroke'
          mark=' m'
        else
          kx=xfac*(xy(j  ) - axes(1))
          ky=yfac*(xy(j+1) - axes(3))
          write (9,*) kx,ky,mark
          mark=' l'
        endif
 1800 continue
      if (mark .eq. ' l') write(9,*)'stroke'
      write(9,*) '0 0 0 sethsbcolor 4 setlinewidth'
      return
c
c  Error conditions
 3100 write(*,'(a,a)')
     $'>>> Error reading x-y line data in file: ',list
      return
      end
c______________________________________________________________________
*====================================
*====================================
* Color palette data base BEGIN=================
*
      blockdata artist
c  Sets the color attributes of the various palette and default color
c  schemes.  Setting the variables in /pallet/ is the main task
c
      common /local/ color(1),black(2)
c
      parameter (levmx=100)
      common /noir/ scutum,kolors,art(4,levmx),kart(levmx),loutln(levmx)
c
      common /pallet/ many,npal(10),hue(28,10),sat(28,10),bri(28,10)
c
      data color/0/,  black/ 0.7, 2.0/
      data kolors/0/, art/levmx*1.0,levmx*1.0,levmx*1.0,levmx*1.0/
c
      data many/7/
c
c  Palette 1:  Blue-Red.  low cool blue to high hot red
      data npal(1)/16 /
      data (hue(j,1),sat(j,1),bri(j,1), j=1, 16)/
     $0.600,1.000,0.500, 0.600,1.000,0.600, 0.600,1.000,0.700,
     $0.600,1.000,1.000, 0.600,0.800,1.000, 0.600,0.700,1.000,
     $0.600,0.500,1.000, 0.600,0.300,1.000, 1.000,0.100,1.000,
     $1.000,0.300,1.000, 1.000,0.350,1.000, 0.900,0.500,1.000,
     $1.000,0.550,1.000, 1.000,0.600,1.000, 1.000,1.000,1.000,
     $1.000,0.500,0.800/
c
c  Palette 2: Hysteria: bright  blue-green-yellow-orange-red-brown
      data npal(2)/24 /
      data (hue(j,2),sat(j,2),bri(j,2), j=1, 24)/
     $0.75,0,0,
     $0.750,1.000,0.5,
     $0.750,1.000,0.7,
     $0.700,1.000,1.000,
     $0.550,1.000,0.800,
     $0.500,1.000,0.700,
     $0.450,1.000,0.800,
     $0.400,1.000,0.900,
     $0.300,0.600,1.000,
     $0.250,0.600,1.000,
     $0.200,0.600,1.000,
     $0.150,0.300,1.000,
     $0.100,0.400,1.000,
     $0.100,0.500,1.000,
     $0.100,0.600,1.000,
     $0.100,0.800,1.000,
     $0.100,1.000,1.000,
     $0.000,0.500,1.000,
     $0.000,0.500,0.900,
     $0.000,1.000,1.000,
     $0.000,1.000,0.800,
     $0.000,1.000,0.600,
     $0.000,1.000,0.500,
     $0.000,0.500,0.300/
c
c  Palette 3: Rainbow = Full sat blue-green-yellow-orange-red
      data npal(3)/7/
      data (hue(j,3),sat(j,3),bri(j,3), j=1, 7)/
     $ 0.7, 1, 1,
     $ 0.6, 1, 1,
     $ 0.5, 1, 1,
     $ 0.2, 1, 1,
     $ 0.15,1, 1,
     $ 0.1, 1, 1,
     $ 0.0, 1, 1/
c
c  Palette 4: Harvard = dark blue-light blue-yellow-brown
      data npal(4)/13/
      data (hue(j,4),sat(j,4),bri(j,4), j=1, 13)/
     $0.65,0.9,0.5,
     $0.65,0.9,0.7,
     $0.65,0.9,1,
     $0.6,0.8,1,
     $0.55,0.7,1,
     $0.55,0.5,1,
     $0.55,0.3,1,
     $0.10,0.3,1,
     $0.10,0.5,1,
     $0.05,0.5,1,
     $0.05,0.9,1,
     $0.05,0.90,0.7,
     $0.05,0.9,0.5/
c
c  Palette 5: Sand = yellow-orange-lightbrown
      data npal(5)/12/
      data (hue(j,5),sat(j,5),bri(j,5), j=1, 12)/
     $ 0.150,0.300,1.000,
     $ 0.145,0.364,1.000,
     $ 0.141,0.427,1.000,
     $ 0.136,0.491,1.000,
     $ 0.132,0.555,1.000,
     $ 0.127,0.618,1.000,
     $ 0.123,0.682,1.000,
     $ 0.118,0.745,1.000,
     $ 0.114,0.809,1.000,
     $ 0.109,0.873,1.000,
     $ 0.105,0.936,1.000,
     $ 0.100,1.000,1.000/
c
c  Palette 6: Cartographic = pale blue-green-yellow-orange 
      data npal(6)/7/
      data (hue(j,6),sat(j,6),bri(j,6), j=1,  7)/
     $ 0.7, 0.5, 1,
     $ 0.6, 0.5, 1,
     $ 0.5, 0.5, 1,
     $ 0.2, 0.5, 1,
     $ 0.15,0.5, 1,
     $ 0.1, 0.5, 1,
     $ 0.0, 0.5, 1/
c  Palette 7: Oceanographic = Dark blue - light blue
      data npal(7)/ 8/
      data (hue(j,7),sat(j,7),bri(j,7), j=1, 8)/
     $ 0.65, 1, 0.5,
     $ 0.65, 1, 0.7,
     $ 0.65, 1, 1,
     $ 0.60, 1, 1,
     $ 0.55, 1, 1,
     $ 0.50, 1, 1,
     $ 0.50, 0.5, 1,
     $ 0.50, 0.2, 1/
c
      end
c_______________________________________________________________
*====================== Color palette data base END
      subroutine encaps(iact, x1,y1,x2,y2)
c$$$$ calls nothing
c  Keeps track of bounding rectangle to create EPS numbers
c  which are written at the end of the plotfile after showpage
      dimension bbox(4)
      save bbox,xorg,yorg
      data bbox/1000,1000,-1000,-1000/,xorg,yorg/0,0/
c
      goto (1000,2000,3000 ),iact
c
c  Origin has moved
 1000 xorg=x1 + xorg
      yorg=y1 + yorg
      return
c
c  Extend bounding box as necessary
 2000 bbox(1) = min(bbox(1), min(x1,x2)+xorg)
      bbox(2) = min(bbox(2), min(y1,y2)+yorg)
      bbox(3) = max(bbox(3), max(x1,x2)+xorg)
      bbox(4) = max(bbox(4), max(y1,y2)+yorg)
      return
c
c  Write out bounding box at end of plotfile
 3000 write(9,'(a/a,4i8)') '%!PS-Adobe-2.0 EPSF-1.2',
     $   '%%BoundingBox:',(nint(72*bbox(j)),j=1,4)
      return
c
      end
c_______________________________________________________________________
