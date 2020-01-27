*CMZ :          28/11/2003  09.55.49  by  Unknown
*-- Author :    Unknown   28/11/2003


      subroutine dvcs_end


* --- Close the ntuple unit
      CALL NTEND

* --- Write the file dvcs_bases_output.txt
         LU     = 66
         CALL SPINFO( LU )
         CALL SHPLOT( LU )
        close(66)


* --- Close bases.data
      close(23)

      return
      end
