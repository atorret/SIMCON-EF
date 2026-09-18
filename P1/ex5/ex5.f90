PROGRAM EX5
    IMPLICIT NONE
    REAL :: m=0.200,k=2.0
    REAL :: t_max=10.0, dt, x_0=0.05, v_0=0.0, t_0=0.0
    REAL, ALLOCATABLE :: x_euler_001(:),v_euler_001(:), t_euler_001(:)
    REAL, ALLOCATABLE :: x_euler_01(:),v_euler_01(:), t_euler_01(:)
    REAL, ALLOCATABLE :: x_euler_02(:),v_euler_02(:), t_euler_02(:)

    REAL, ALLOCATABLE :: x_euler_predictor(:),v_euler_predictor(:), t_euler_predictor(:)
    REAL, ALLOCATABLE :: x_verlet(:),v_verlet(:), t_verlet(:)

    dt=0.001
    CALL EULER(m,k,t_max,dt,x_0,v_0,t_0,x_euler_001,v_euler_001,t_euler_001)
    dt=0.01
    CALL EULER(m,k,t_max,dt,x_0,v_0,t_0,x_euler_01,v_euler_01,t_euler_01)
    dt=0.02
    CALL EULER(m,k,t_max,dt,x_0,v_0,t_0,x_euler_02,v_euler_02,t_euler_02)
    CALL EULER_PREDICTOR(m,k,t_max,dt,x_0,v_0,t_0,x_euler_predictor,v_euler_predictor,t_euler_predictor)
    CALL VERLET(m,k,t_max,dt,x_0,v_0,t_0,x_verlet,v_verlet,t_verlet)

CONTAINS
    SUBROUTINE EULER(m,k,t_max,dt,x_0,v_0,t_0,x,v,t)
        IMPLICIT NONE
        REAL, INTENT(IN) :: m,k,t_max,dt,x_0,v_0,t_0
        REAL, DIMENSION(:), INTENT(OUT) :: x,v,t
        REAL :: a
        INTEGER :: i

        x(1)=x_0
        v(1)=v_0
        t(1)=t_0

        DO i = 2, INT(t_max/dt)+1
            a=-(k/m)*x(i)
            v(i)=v(i-1)+a*dt
            x(i)=x(i-1)+v(i-1)*dt+0.5*a*dt**2
            t(i)=t(i-1)+dt
        END DO
    END SUBROUTINE EULER

    SUBROUTINE VERLET(m,k,t_max,dt,x_0,v_0,t_0,x,v,t)
        IMPLICIT NONE
        REAL, INTENT(IN) :: m,k,t_max,dt,x_0,v_0,t_0
        REAL, DIMENSION(:), INTENT(OUT) :: x,v,t
        REAL :: a
        INTEGER :: i

        x(1)=x_0
        v(1)=v_0
        t(1)=t_0

        DO i = 2, INT(t_max/dt)+1
            a=-(k/m)*x(i)
            v(i)=v(i-1)+a*dt
            x(i)=x(i-1)+v(i-1)*dt+0.5*a*dt**2
            t(i)=t(i-1)+dt
        END DO
        RETURN
    END SUBROUTINE VERLET

    SUBROutine EULER_PREDICTOR(m,k,t_max,dt,x_0,v_0,t_0,x,v,t)
        IMPLICIT NONE
        REAL, INTENT(IN) :: m,k,t_max,dt,x_0,v_0,t_0
        REAL, DIMENSION(:), INTENT(OUT) :: x,v,t
        REAL :: a
        INTEGER :: i

        x(1)=x_0
        v(1)=v_0
        t(1)=t_0

        DO i = 2, INT(t_max/dt)+1
            a=-(k/m)*x(i)
            v(i)=v(i-1)+a*dt
            x(i)=x(i-1)+v(i-1)*dt+0.5*a*dt**2
            t(i)=t(i-1)+dt
        END DO
        RETURN
    END SUBROUTINE EULER_PREDICTOR

END PROGRAM EX5