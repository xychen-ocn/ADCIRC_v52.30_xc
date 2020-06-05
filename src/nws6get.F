C***********************************************************************
C                                                                      *
C                                                                      *
C                                                                      *
C   Read in and interpolate onto the ADCIRC grid wind and pressure     *
C   fields from a meteorological file on a rectangular grid (either in *
C   Longitude, Latitude or Cartesian coordinates, consistent with the  *
C   ADCIRC grid coordinates).  If the ADCIRC grid is in Lon,Lat these  *
C   MUST BE IN RADIANS!                                                *
C                                                                      *
C   It is assumed that the meteorological grid is set up so that y     *
C   (e.g., latitude) varies from north (k=1) to south (k=NWLAT) and x  *
C   (e.g., longitude) varies from west (j=1) to east (j=NWLON).        *
C                                                                      *
C   The spatial extents of the meteorological grid must be consistent  *
C   with the ADCIRC model domain.  For example, if ADCIRC uses negative*
C   longitude values to indicate locations W of the Greenwich meridian,*
C   the meteorological file must be similarly organized.  Any grid that*
C   crosses the Greenwich Meridian should be organized so that the seam*
C   occurs @ 180 deg longitude. Therefore, the meteorological and      *
C   ADCIRC grids should use negative longitudes W of the Greenwich     *
C   Meridian and positive longitudes to the E.                         *
C                                                                      *
C                                                                      *
C   NOTE:  It is assumed that the met file data is oriented so that    *
C          the outer loop is on latitude and the inner loop is on      *
C          longitude.  For example:                                    *
C          line 1             lat 1,     lon 1                         *
C          line 2             lat 1,     lon 2                         *
C            .                                                         *
C          line nwlon         lat 1,     lon nwlon                     *
C          line nwlon+1       lat 2,     lon 1                         *
C          line nwlon+2       lat 2,     lon 2                         *
C            .                                                         *
C          line 2*nwlon       lat 2,     lon nwlon                     *
C          line 2*nwlon+1     lat 3,     lon 1                         *
C          line 2*nwlon+2     lat 3,     lon 2                         *
C            .                                                         *
C          line nwlon*nwlat   lat nwlat, lon nwlon                     *
C                                                                      *
C   NOTE:  It is assumed that he met file data is oriented so that     *
C          latitude varies from the northern most value (lat 1) to the *
C          southern most value (lat nwlat) and longitude varies in an  *
C          easterly direction (e.g. from 0 to 360 where positive       *
C          longitudes are angles measured easterly of the GM.          *
C                                                                      *
C   NOTE:  For the global AVN grid running from 0.5 - 359.5 deg lon    *
C          and 90 - -90 deg lat in 1 degree increments, NWLAT=181 and  *
C          NWLON=360 yielding a total number of entries in the file    *
C          of 65160.                                                   *
C                                                                      *
C   NOTE:  It is assumed that wind velocity is in EAST,NORTH components*
C          in M/2 and pressure is in N/M^2                             *
C                                                                      *
C   NOTE:  WLATMAX,WLONMIN,WLATINC,WLONINC should be in deg.           *
C                                                                      *
C   NOTE:  This should wrap if XCOORD > WLONMIN+NWLON*WLONINC  or      *
C          XCOORD < WLONMIN                                            *
C                                                                      *
C                                                                      *
C   MNWLAT = MAXIMUM NUMBER OF LATITUDES IN WIND FILE                  *
C            SET = 1 IF FLEET NUMERIC WIND FILE NOT IN USE             *
C   MNWLON = MAXIMUM NUMBER OF LONGITUDES IN WIND FILE                 *
C            SET = 1 IF FLEET NUMERIC WIND FILE NOT IN USE             *
C                                                                      *
C                           R.L. 4/13/99                               *
C                           R.L.09/04/00 added RHOWAT0 to call         *
C   R.L.09/04/00 added RHOWAT0 to call                                 *
C   R.L. 4/2/01  changed MNWLAT,MNWLON in ALLOCATE statement to        *
C                NWLAT,NWLON                                           *
C   R.L. 8/10/05 eliminated adding 360 to negative longitudes to match *
C                AVN model grid setup.  User is now required to provide*
C                met and ADCIRC grid that correspond in space without  *
C                adjusting the longitude.  Also the input variable     *
C                order has been changed to U,V,P to be consistent with *
C                other NWS input formats.                              *
C***********************************************************************

      SUBROUTINE NWS6GET(X,Y,SLAM,SFEA,WVNX,WVNY,PRESS,NP,NWLON,NWLAT,
     &     WLATMAX,WLONMIN,WLATINC,WLONINC,ICS,RHOWAT0,G)
      USE SIZES

      IMPLICIT NONE
      INTEGER, SAVE :: FIRSTCALL = 0
      INTEGER NP,NWLON,NWLAT,I,J,ICS
      REAL(SZ) RHOWAT0,RHOWATG,G
      INTEGER  LATIND1,LATIND2,LONIND1,LONIND2
      REAL(SZ) WLATMAX,WLONMIN,WLATINC,WLONINC,XWRATIO,YWRATIO
      REAL(SZ) WLATM,WLONM
      REAL*8 X(*),Y(*),SLAM(*),SFEA(*),XCOOR,YCOOR
      REAL(SZ) WVNX(*),WVNY(*),PRESS(*)
      REAL(SZ),SAVE,ALLOCATABLE :: WVXFN(:,:),WVYFN(:,:),PRN(:,:)
C
      IF (FIRSTCALL.EQ.0) THEN
         FIRSTCALL = 1
         ALLOCATE ( WVXFN(NWLAT,NWLON),WVYFN(NWLAT,NWLON),
     &        PRN(NWLAT,NWLON) )
      ENDIF
C
      RHOWATG=RHOWAT0*G
      DO I=1,NWLAT
         DO J=1,NWLON
            READ(22,*) WVXFN(I,J),WVYFN(I,J),PRN(I,J)
         END DO
      END DO

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
      END SUBROUTINE NWS6GET
