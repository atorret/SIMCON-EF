PROGRAM EX5
    IMPLICIT NONE
    REAL :: m = 0.200, k = 2
    REAL :: t_max = 10.0, x_0 = 0.05, v_0 = 0.0, t_0 = 0.0

    CALL EULER(m, k, t_max, 0.001, x_0, v_0, t_0, "euler_001.dat")
    CALL EULER(m, k, t_max, 0.010, x_0, v_0, t_0, "euler_01.dat")
    CALL EULER(m, k, t_max, 0.020, x_0, v_0, t_0, "euler_02.dat")
    CALL EULER_PREDICTOR(m, k, t_max, 0.020, x_0, v_0, t_0, "euler_predictor.dat")
    CALL VERLET(m, k, t_max, 0.020, x_0, v_0, t_0, "verlet.dat")

    PRINT *, "Wrote euler_*.dat, euler_predictor.dat, verlet.dat"

CONTAINS

    SUBROUTINE EULER(m, k, t_max, dt, x_0, v_0, t_0, fname)
        IMPLICIT NONE
        REAL, INTENT(IN) :: m, k, t_max, dt, x_0, v_0, t_0
        CHARACTER(LEN=*), INTENT(IN) :: fname
        REAL :: x, v, t, a, ek, ep, et
        INTEGER :: i, n, u

        n = INT((t_max - t_0) / dt)
        x = x_0
        v = v_0
        t = t_0
        u = 20

        OPEN(UNIT=u, FILE=fname, STATUS="REPLACE", ACTION="WRITE")
        WRITE(u, '(A)') "# t x v Ek Ep ET"
        WRITE(u, '(6ES16.8)') t, x, v, 0.5*m*v**2, 0.5*k*x**2, 0.5*(m*v**2 + k*x**2)

        DO i = 1, n
            a = -(k/m) * x
            x = x + v*dt + 0.5*a*dt**2
            v = v + a*dt
            t = t + dt
            ek = 0.5*m*v**2
            ep = 0.5*k*x**2
            et = ek + ep
            WRITE(u, '(6ES16.8)') t, x, v, ek, ep, et
        END DO
        CLOSE(u)
    END SUBROUTINE EULER

    SUBROUTINE VERLET(m, k, t_max, dt, x_0, v_0, t_0, fname)
        IMPLICIT NONE
        REAL, INTENT(IN) :: m, k, t_max, dt, x_0, v_0, t_0
        CHARACTER(LEN=*), INTENT(IN) :: fname
        REAL :: x, v, t, a, a_new, ek, ep, et
        INTEGER :: i, n, u

        n = INT((t_max - t_0) / dt)
        x = x_0
        v = v_0
        t = t_0
        a = -(k/m) * x
        u = 20

        OPEN(UNIT=u, FILE=fname, STATUS="REPLACE", ACTION="WRITE")
        WRITE(u, '(A)') "# t x v Ek Ep ET"
        WRITE(u, '(6ES16.8)') t, x, v, 0.5*m*v**2, 0.5*k*x**2, 0.5*(m*v**2 + k*x**2)

        DO i = 1, n
            x = x + v*dt + 0.5*a*dt**2
            a_new = -(k/m) * x
            v = v + 0.5*(a + a_new)*dt
            a = a_new
            t = t + dt
            ek = 0.5*m*v**2
            ep = 0.5*k*x**2
            et = ek + ep
            WRITE(u, '(6ES16.8)') t, x, v, ek, ep, et
        END DO
        CLOSE(u)
    END SUBROUTINE VERLET

    SUBROUTINE EULER_PREDICTOR(m, k, t_max, dt, x_0, v_0, t_0, fname)
        IMPLICIT NONE
        REAL, INTENT(IN) :: m, k, t_max, dt, x_0, v_0, t_0
        CHARACTER(LEN=*), INTENT(IN) :: fname
        REAL :: x, v, t, a, x_p, a_p, ek, ep, et
        INTEGER :: i, n, u

        n = INT((t_max - t_0) / dt)
        x = x_0
        v = v_0
        t = t_0
        u = 20

        OPEN(UNIT=u, FILE=fname, STATUS="REPLACE", ACTION="WRITE")
        WRITE(u, '(A)') "# t x v Ek Ep ET"
        WRITE(u, '(6ES16.8)') t, x, v, 0.5*m*v**2, 0.5*k*x**2, 0.5*(m*v**2 + k*x**2)

        DO i = 1, n
            a = -(k/m) * x
            x_p = x + v*dt + 0.5*a*dt**2
            a_p = -(k/m) * x_p
            a = 0.5*(a_p + a)
            x = x + v*dt + 0.5*a*dt**2
            v = v + a*dt
            t = t + dt
            ek = 0.5*m*v**2
            ep = 0.5*k*x**2
            et = ek + ep
            WRITE(u, '(6ES16.8)') t, x, v, ek, ep, et
        END DO
        CLOSE(u)
    END SUBROUTINE EULER_PREDICTOR

END PROGRAM EX5