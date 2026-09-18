PROGRAM EX5
    IMPLICIT NONE
    REAL :: m,k,x,v,x_o,v_o
    REAL :: t, t_max, dt
    m=0.200
    k=2.0
    x_o=0.05
    v_o=0.0

    t=0.0
    t_max=10-0

    dt=0.001
    CALL EULER(m,k,x,v,x_o,v_o,t,t_max,dt)
    dt=0.01
    CALL EULER(m,k,x,v,x_o,v_o,t,t_max,dt)
    dt=0.02
    CALL EULER(m,k,x,v,x_o,v_o,t,t_max,dt)
    CALL EULER_PREDICTOR(m,k,x,v,x_o,v_o,t,t_max,dt)
    CALL VERLET(m,k,x,v,x_o,v_o,t,t_max,dt)

SUBROUTINE EULER(m,k,x,v,x_o,v_o,t,t_max,dt)
    IMPLICIT NONE
    REAL :: m,k,x,v,x_o,v_o
    REAL :: t, t_max, dt
    REAL :: a
    a=-(k/m)*x
    v=v_o+a*dt
    x=x_o+v*dt
    RETURN
END SUBROUTINE EULER

SUBROUTINE VERLET(m,k,x,v,x_o,v_o,t,t_max,dt)
    IMPLICIT NONE
    REAL :: m,k,x,v,x_o,v_o
    REAL :: t, t_max, dt
    REAL :: a
    a=-(k/m)*x
    v=v_o+a*dt
    x=x_o+v*dt
    RETURN
END SUBROUTINE VERLET

SUBROutine EULER_PREDICTOR(m,k,x,v,x_o,v_o,t,t_max,dt)
    IMPLICIT NONE
    REAL :: m,k,x,v,x_o,v_o
    REAL :: t, t_max, dt
    REAL :: a
    a=-(k/m)*x
    v=v_o+a*dt
    x=x_o+v*dt
    RETURN
END SUBROUTINE EULER_PREDICTOR


END PROGRAM EX5