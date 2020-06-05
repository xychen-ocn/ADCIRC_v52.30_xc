C     ----------------------------------------------------------------
C      S U B R O U T I N E     N W S 7 G E T
C     ----------------------------------------------------------------
C
C     jgf46.01 Subroutine to get a set of surface wind stresses and
C     barometric pressure on a rectangular grid (either in Longitude,
C     Latitude or Cartesian coordinates, consistent with the ADCIRC grid
C     coordinates) from unit 22 and interpolate them in space onto the
C     ADCIRC grid. If the ADCIRC grid is in Lon, Lat these must be in
C     radians.
C
C     It is assumed that the meteorological grid is set up so that y
C     (e.g., latitude) varies from north (k=1) to south (k=NWLAT) and x
C     (e.g., longitude) varies from west (j=1) to east (j=NWLON).
C
C     The spatial extents of the meteorological grid must be consistent
C     with the ADCIRC model domain.  For example, if ADCIRC uses
C     negative longitude values to indicate locations W of the Greenwich
C     meridian, the meteorological file must be similarly organized.
C     Any grid that crosses the Greenwich Meridian should be organized
C     so that the seam occurs @ 180 deg longitude. Therefore, the
C     meteorological and ADCIRC grids should use negative longitudes W
C     of the Greenwich Meridian and positive longitudes to the E.
C
C   NOTE:  It is assumed that the met file data is oriented so that
C          the outer loop is on latitude and the inner loop is on
C          longitude.  For example:
C          line 1             lat 1,     lon 1
C          line 2             lat 1,     lon 2
C            .
C          line nwlon         lat 1,     lon nwlon
C          line nwlon+1       lat 2,     lon 1
C          line nwlon+2       lat 2,     lon 2
C            .
C          line 2*nwlon       lat 2,     lon nwlon
C          line 2*nwlon+1     lat 3,     lon 1
C          line 2*nwlon+2     lat 3,     lon 2
C            .
C          line nwlon*nwlat   lat nwlat, lon nwlon
C
C   NOTE:  It is assumed that he met file data is oriented so that
C          latitude varies from the northern most value (lat 1) to the
C          southern most value (lat nwlat) and longitude varies in an
C          easterly direction (e.g. from 0 to 360 where positive
C          longitudes are angles measured easterly of the GM.
C
C   NOTE:  It is assumed that wind stress is in EAST, NORTH components
C          in m/s and pressure is in N/m^2
C
C   NOTE:  WLATMAX,WLONMIN,WLATINC,WLONINC should be in deg.
C
C   NOTE:  This should wrap if XCOORD > WLONMIN+NWLON*WLONINC  or
C          XCOORD < WLONMIN
C
C
C
C
C
C
C
C
C
C
C
C
C
C
C
C
C
C
C
C
C
C
C
C     ----------------------------------------------------------------

      SUBROUTINE NWS7GET(X,Y,SLAM,SFEA,WVNX,WVNY,PRESS,NP,NWLON,NWLAT,
     &     WLATMAX,WLONMIN,WLATINC,WLONINC,ICS,RHOWAT0,G)
      IMPLICIT NONE
      INTEGER, intent(in) :: NP     ! number of nodes
      REAL(8), intent(in), dimension(NP) :: X
      REAL(8), intent(in), dimension(NP) :: Y
      REAL(8), intent(in), dimension(NP) :: SLAM
      REAL(8), intent(in), dimension(NP) :: SFEA
      REAL(SZ), intent(out), dimension(NP) :: WVNX
      REAL(SZ), intent(out), dimension(NP) :: WVNY
      REAL(SZ), intent(out), dimension(NP) :: PRESS
      INTEGER, intent(in) :: NWLON  ! max. # of longitudes in stress file
      INTEGER, intent(in) :: NWLAT  ! max. # of latitutes in stress file
      REAL(SZ), intent(in) :: WLATMAX
      REAL(SZ), intent(in) :: WLONMIN
      REAL(SZ), intent(in) :: WLATINC
      REAL(SZ), intent(in) :: WLONINC
      INTEGER, intent(in) :: ICS    ! coord.sys., 1=cartesian, 2=spherical
      REAL(SZ), intent(in):: RHOWAT0! ref. dens. of water
      REAL(SZ), intent(in):: G      ! gravitational constant

      INTEGER I                     ! node loop counter
      INTEGER K, J                  ! latitude, longitude loop counters
      REAL(SZ) RHOWATG             ! ref. dens. of water * grav. constant
      INTEGER  LATIND1,LATIND2,LONIND1,LONIND2
      REAL(SZ) WLATM,WLONM,XWRATIO,YWRATIO
      REAL(8) XCOOR,YCOOR
      REAL(SZ), SAVE, ALLOCATABLE :: WVXFN(:,:),WVYFN(:,:),PRN(:,:)
      LOGICAL, SAVE :: MemoryAllocated = .False.
C
      IF (.not.MemoryAllocated) THEN
         ALLOCATE ( WVXFN(NWLAT,NWLON),WVYFN(NWLAT,NWLON),
     &        PRN(NWLAT,NWLON) )
         MemoryAllocated = .True.
      ENDIF
C
      RHOWATG=RHOWAT0*G
      DO K=1,NWLAT
         DO J=1,NWLON
            READ(22,*) WVXFN(K,J),WVYFN(K,J),PRN(K,J)
         END DO
      END DO
C
      DO I=1,NP                 !INTERPOLATE TO ADCIRC GRID
         IF(ICS.EQ.2) THEN
            YCOOR=SFEA(I)*RAD2DEG
            XCOOR=SLAM(I)*RAD2DEG
         ENDIF
         IF(ICS.EQ.1) THEN
            YCOOR=Y(I)
            XCOOR=X(I)
         ENDIF
         LATIND2=(WLATMAX-YCOOR)/WLATINC + 1
         IF(LATIND2.EQ.NWLAT) LATIND2=LATIND2-1
         LATIND1=LATIND2 + 1
         LONIND1=(XCOOR-WLONMIN)/WLONINC + 1
         LONIND2=LONIND1 + 1
C
         WLONM = WLONMIN + (LONIND1-1)*WLONINC
         WLATM = WLATMAX - (LATIND1-1)*WLATINC
         XWRATIO=(XCOOR-WLONM)/WLONINC
         YWRATIO=(YCOOR-WLATM)/WLATINC
C
         IF(LONIND1.EQ.0) LONIND1=NWLON
         IF(LONIND1.EQ.NWLON) LONIND2=1
C
         WVNX(I) = WVXFN(LATIND2,LONIND2)*XWRATIO*YWRATIO
     &        + WVXFN(LATIND2,LONIND1)*(1.d0-XWRATIO)*YWRATIO
     &        + WVXFN(LATIND1,LONIND2)*XWRATIO*(1.d0-YWRATIO)
     &        + WVXFN(LATIND1,LONIND1)*(1.d0-XWRATIO)*(1.d0-YWRATIO)
         WVNY(I) = WVYFN(LATIND2,LONIND2)*XWRATIO*YWRATIO
     &        + WVYFN(LATIND2,LONIND1)*(1.d0-XWRATIO)*YWRATIO
     &        + WVYFN(LATIND1,LONIND2)*XWRATIO*(1.d0-YWRATIO)
     &        + WVYFN(LATIND1,LONIND1)*(1.d0-XWRATIO)*(1.d0-YWRATIO)
         PRESS(I) = PRN(LATIND2,LONIND2)*XWRATIO*YWRATIO
     &        + PRN(LATIND2,LONIND1)*(1.d0-XWRATIO)*YWRATIO
     &        + PRN(LATIND1,LONIND2)*XWRATIO*(1.d0-YWRATIO)
     &        + PRN(LATIND1,LONIND1)*(1.d0-XWRATIO)*(1.d0-YWRATIO)
         PRESS(I) = PRESS(I)/RHOWATG
      END DO
C
      RETURN
C     ----------------------------------------------------------------
      END SUBROUTINE NWS7GET
C     ----------------------------------------------------------------
