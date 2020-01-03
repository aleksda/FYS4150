from planetes import *
import time

class System:

    def __init__(self, planets = None):

        self.planets = []
        self.x0 = []
        self.v0 = []
        self.m = []

        if planets is not None:

            self.__add__(planets)

    # The __add__ function checks if argument planets are lone planet or list/tuple of planets then adds them to the system
    def __add__(self, planets):

        if isinstance(planets, (list, tuple)):
            for p in planets:
                self.planets.append(p)
                self.x0.append(p.x0)
                self.v0.append(p.v0)
                self.m.append(p.m)

        else:
            self.planets.append(planets)
            self.x0.append(planets.x0)
            self.v0.append(planets.v0)
            self.m.append(planets.m)

    def __iadd__(self, planets):
        self.__add__(planets)

    #Solve function that also contains function for solving acceleration(a_g), forward_euler, euler_cromer and velocity_verlet
    def solve(self, T, dt, origin_mass, method = "forward_euler"):

        #Constants that will be used in functions a_g, forward_euler, euler_cromer and velocity_verlet
        G = 4*np.pi**2

        steps = int(T/dt)
        if steps*dt < T:
            steps += 1
        T = steps*dt

        x0 = np.array(self.x0)
        v0 = np.array(self.v0)

        N = x0.shape[0]
        p = x0.shape[1]

        t = np.arange(0, T+dt, dt)
        x = np.zeros((steps, N, p))
        v = np.zeros((steps, N, p))
        m = np.zeros(N)[:,np.newaxis]

        x[0] = x0
        v[0] = v0

        #Function for calculating acceleration. Used by forward_euler, euler_cromer and velocity_verlet
        @njit
        def a_g(m, r, r_norm, G):
            idx = r_norm != 0
            a = G*m[idx]*r[idx]/r_norm[idx].reshape(-1,1)**3
            return np.sum(a, axis = 0)


        #All three functions, have built in procentcounters which prints loading process in terminal
        @njit
        def forward_euler(steps, N, m, x, v, dt, G, origin_mass):
            total = (steps-1)*N
            count = 0
            perc = 0
            perc_new = 0

            for i in range(steps-1):
                for j in range(N):
                    r = x[i] - x[i,j]
                    r_norm = np.sqrt(np.sum(r**2, axis = 1))
                    v[i+1,j] = v[i,j] + a_g(m, r, r_norm, G)*dt
                    if origin_mass != 0:
                        x_norm = np.sqrt(np.sum(x[i,j]**2, axis = 0))
                        v[i+1,j] = v[i+1,j] - G*origin_mass*x[i,j]*dt/x_norm**3
                    x[i+1,j] = x[i,j] + v[i,j]*dt

                    count += 1
                    new_perc = int(100*count/total)
                    if new_perc > perc:
                        perc = new_perc
                        print(perc)
            return x,v

        @njit
        def euler_cromer(steps, N, m, x, v, dt, G, origin_mass):
            total = (steps-1)*N
            count = 0
            perc = 0
            perc_new = 0

            for i in range(steps-1):
                for j in range(N):
                    r = x[i] - x[i,j]
                    r_norm = np.sqrt(np.sum(r**2, axis = 1))
                    v[i+1,j] = v[i,j] + a_g(m, r, r_norm, G)*dt
                    if origin_mass != 0:
                        x_norm = np.sqrt(np.sum(x[i,j]**2, axis = 0))
                        v[i+1,j] = v[i+1,j] - G*origin_mass*x[i,j]*dt/x_norm**3
                    x[i+1,j] = x[i,j] + v[i+1,j]*dt

                    count += 1
                    new_perc = int(100*count/total)
                    if new_perc > perc:
                        perc = new_perc
                        print(perc)
            return x,v

        @njit(fastmath=True, parallel=False)
        def velocity_verlet(steps, N, m, x, v, dt, G, origin_mass):
            total = (steps-1)*N
            count = 0
            perc = 0
            perc_new = 0

            for i in prange(steps-1):
                for j in prange(N):
                    r = x[i] - x[i,j]
                    r_norm = np.sqrt(np.sum(r**2, axis = 1))
                    v_half = v[i,j] + 0.5*a_g(m, r, r_norm, G)*dt
                    if origin_mass != 0:
                        x_norm = np.sqrt(np.sum(x[i,j]**2, axis = 0))
                        v_half = v_half - G*origin_mass*x[i,j]*dt/x_norm**3
                    x[i+1,j] = x[i,j] + v_half*dt

                    r = x[i+1] - x[i+1,j]
                    r_norm = np.sqrt(np.sum(r**2, axis = 1))
                    v[i+1,j] = v_half + a_g(m, r, r_norm, G)*dt
                    if origin_mass != 0:
                        x_norm = np.sqrt(np.sum(x[i+1,j]**2, axis = 0))
                        v[i+1,j] = v[i+1,j] - G*origin_mass*x[i+1,j]*dt/x_norm**3

                    count += 1
                    new_perc = int(100*count/total)
                    if new_perc > perc:
                        perc = new_perc
                        print(perc)
            return x,v

        if method == "forward_euler":
            x,v = forward_euler(steps, N, m, x, v, dt, G, origin_mass)
        elif method == "euler_cromer":
            x,v = euler_cromer(steps, N, m, x, v, dt, G, origin_mass)
        elif method == "velocity_verlet":
            x,v = velocity_verlet(steps, N, m, x, v, dt, G, origin_mass)

        self.t = t
        self.x = x
        self.v = v

    #Function for plotting 3D plots
    def plot(self, max_points = 1E4):
        steps = self.x.shape[0]
        skip = np.min((np.max((int(steps//max_points), 1)), steps))
        fig = plt.figure(figsize=plt.figaspect(1)*1.5)
        ax = fig.gca(projection='3d')
        ax.scatter(xs=0, ys=0, zs=0, c='C0', label="Sun")
        ax.set_xlabel('x [AU]')
        ax.set_ylabel('y [AU]')
        ax.set_zlabel('z [AU]')


        d_max = np.max([np.max(np.abs(self.x[:,:,0])),
                            np.max(np.abs(self.x[:,:,1])),
                            np.max(np.abs(self.x[:,:,2]))])

        for i in range(len(self.planets)):
            ax.plot(self.x[::skip,i,0], self.x[::skip,i,1], self.x[::skip,i,2],
                    label = self.planets[i].name, c=f'C{i+1}')
        ax.legend()
        ax.set_xlim([-d_max, d_max])
        ax.set_ylim([-d_max, d_max])
        ax.set_zlim([-d_max, d_max])
        plt.show()

if __name__ == "__main__":
    
    t0 = time.process_time()

    # oppg. C
    # to make circular orbits
    #se = System(Mercury)
    se = System(Earth)
    se.solve(1, 1E-4, 1, method = "velocity_verlet")
    #se.solve(1, 1E-4, 1, method = "forward_euler")
    t1 = time.process_time()
    print(f"Calculations took {t1-t0} sec")
    se.plot()
    
    """
    # to find breakpoints as function of dt
    se.solve(0.5, 0.01, 1, method = "velocity_verlet")
    #se.solve(2, 0.01, 1, method = "forward_euler")
    se.plot()
    """

    
    # oppg. E
    """
    three = [Earth, Jupiter]
    t = System(three)
    t.solve(40, 1E-4, 1, method="velocity_verlet")
    t.plot()
    """
    """
    # oppg. F
    planets = [Mercury, Venus, Earth, Mars, Jupiter, Saturn, Uranus, Neptune, Pluto]
    S = System(planets)
    S.solve(100, 1E-4, 1, method = "velocity_verlet")
    #S.solve(100, 0.01, 1, method = "forward_euler")
    S.plot()
    """