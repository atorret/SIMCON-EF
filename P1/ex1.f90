PROGRAM EX1
    IMPLICIT NONE
    INTEGER :: A, B
    INCLUDE '../chdir_to_code.inc'

    PRINT *, 'Enter 2 numbers:'
    READ(*, *) A, B

    OPEN(UNIT=10, FILE='ex1_results.txt', STATUS='REPLACE', ACTION='WRITE')
    WRITE(10, '(A, I0, A, I0)') 'The integer numbers are: ', A, ', ', B
    CLOSE(10)

END PROGRAM EX1
