PROGRAM EX4
    IMPLICIT NONE
    INTEGER, PARAMETER :: nmax = 1000
    INTEGER, PARAMETER :: nruns = 3
    REAL, DIMENSION(nmax) :: x, y
    REAL :: A, A_average, A_exact, A_deviation, A_deviation_average
    INTEGER :: N, J

    A_exact = 4.0 * (1.0 + LOG(11.0))

    OPEN(UNIT=10, FILE='ex4_results.dat', STATUS='REPLACE', ACTION='WRITE')
    DO N = 100, nmax,  100
        A_average = 0.0
        A_deviation = 0.0
        DO J = 1, nruns
            CALL RANDOM_NUMBER(x(1:N))
            CALL RANDOM_NUMBER(y(1:N))
            x(1:N) = x(1:N) * 2.0 - 1.0
            y(1:N) = y(1:N) * 22.0 - 11.0
            CALL calculate_area(x, y, N, A)
            A_average = A_average + A
            A_deviation = A_deviation + (A - A_exact)**2
        END DO
        A_average = A_average / REAL(nruns)
        A_deviation_average = SQRT(A_deviation / REAL(nruns))
        WRITE(10, '(I8, 1X, ES16.8, 1X, ES16.8)') N, A_average, A_deviation_average
    END DO
    CLOSE(10)

    PRINT *, 'Results written to ex4_results.dat'
    PRINT *, 'Exact area = ', A_exact

CONTAINS

    SUBROUTINE calculate_area(x, y, N, A)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: N
        REAL, DIMENSION(:), INTENT(IN) :: x, y
        REAL, INTENT(OUT) :: A
        INTEGER :: k, Nh

        Nh = 0
        DO k = 1, N
            IF (x(k) /= 0.0) THEN
                IF (ABS(y(k)) <= ABS(1.0 / x(k))) THEN
                    Nh = Nh + 1
                END IF
            END IF
        END DO
        A = 44.0 * REAL(Nh) / REAL(N)
    END SUBROUTINE calculate_area

END PROGRAM EX4
