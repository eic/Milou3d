*CMZ :          28/11/2003  09.55.49  by  Unknown
*-- Author :    Unknown   28/11/2003
c modified June 2011 S.F.
c last seed produced by drn.f saved in fort.77 ascan 
c be used as input for the new
c generation of spring
c not understood common RANDM in bsread renamed RANDM1 to avoid overwrite 
c of the cards seed 

      subroutine dvcs_end

      Integer IDUMY
      COMMON /MCSEED/IDUMY

* --- Close the ntuple unit
      CALL NTEND

* --- Write the file dvcs_bases_output.txt
         LU     = 66
         CALL SPINFO( LU )
         CALL SHPLOT( LU )
        close(66)

        write(77,*) IDUMY

* --- Close bases.data
      close(23)

      return
      end
