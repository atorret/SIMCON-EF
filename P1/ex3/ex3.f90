PROGRAM EX3
    IMPLICIT NONE
    INTEGER, PARAMETER :: npts = 2001
    INTEGER, PARAMETER :: wmax = 500
    REAL, DIMENSION(npts) :: C, t
    REAL, DIMENSION(0:wmax) :: result
    INTEGER :: i, w
    INCLUDE '../../chdir_to_code.inc'

    OPEN(UNIT=10, FILE='ex3_data.dat', STATUS='OLD', ACTION='READ')
    DO i = 1, npts
        READ(10, *) t(i), C(i)
    END DO
    CLOSE(10)

    CALL simpson(C, t, result)

    OPEN(UNIT=20, FILE='ex3_results.dat', STATUS='REPLACE', ACTION='WRITE')
    DO w = 0, wmax
        WRITE(20, '(I8, 1X, ES16.8)') w, result(w)
    END DO
    CLOSE(20)

    PRINT *, 'Results written to ex3_results.dat'
    PRINT *, 'F(w=0) = ', result(0)

CONTAINS

    SUBROUTINE simpson(C, t, result)
        IMPLICIT NONE
        REAL, DIMENSION(:), INTENT(IN) :: C, t
        REAL, DIMENSION(0:), INTENT(OUT) :: result
        INTEGER :: i, w, n
        REAL :: h, s
        REAL, ALLOCATABLE :: f(:)

        n = SIZE(t)
        h = t(2) - t(1)
        ALLOCATE(f(n))

        DO w = 0, 500
            DO i = 1, n
                f(i) = C(i) * COS(REAL(w) * t(i))
            END DO

            ! Composite Simpson (even number of intervals, step 2)
            s = 0.0
            DO i = 2, n - 1, 2
                s = s + h / 3.0 * (f(i-1) + 4.0 * f(i) + f(i+1))
            END DO
            result(w) = s
        END DO

        DEALLOCATE(f)
    END SUBROUTINE simpson

END PROGRAM EX3
