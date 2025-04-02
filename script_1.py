import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation


def quaternion_to_rotation_matrix(q):
    q0, q1, q2, q3 = q
    return np.array([
        [1 - 2 * (q2 ** 2 + q3 ** 2), 2 * (q1 * q2 - q0 * q3), 2 * (q1 * q3 + q0 * q2)],
        [2 * (q1 * q2 + q0 * q3), 1 - 2 * (q1 ** 2 + q3 ** 2), 2 * (q2 * q3 - q0 * q1)],
        [2 * (q1 * q3 - q0 * q2), 2 * (q2 * q3 + q0 * q1), 1 - 2 * (q1 ** 2 + q2 ** 2)]
    ])


def quaternion_multiply(q, p):
    q0, q1, q2, q3 = q
    p0, p1, p2, p3 = p
    return np.array([
        q0 * p0 - q1 * p1 - q2 * p2 - q3 * p3,
        q0 * p1 + q1 * p0 + q2 * p3 - q3 * p2,
        q0 * p2 - q1 * p3 + q2 * p0 + q3 * p1,
        q0 * p3 + q1 * p2 - q2 * p1 + q3 * p0
    ])


def rk4_step(q, w, dt, torque, I, I_inv):
    def euler_eq(w):
        return np.dot(I_inv, torque - np.cross(w, np.dot(I, w)))

    def quaternion_eq(q, w):
        omega = np.array([0, *w])  # Угловая скорость в виде кватерниона
        return 0.5 * quaternion_multiply(q, omega)

    k1_w = euler_eq(w)
    k2_w = euler_eq(w + 0.5 * dt * k1_w)
    k3_w = euler_eq(w + 0.5 * dt * k2_w)
    k4_w = euler_eq(w + dt * k3_w)
    w_new = w + (dt / 6) * (k1_w + 2 * k2_w + 2 * k3_w + k4_w)

    k1_q = quaternion_eq(q, w)
    k2_q = quaternion_eq(q + 0.5 * dt * k1_q, w)
    k3_q = quaternion_eq(q + 0.5 * dt * k2_q, w)
    k4_q = quaternion_eq(q + dt * k3_q, w)
    q_new = q + (dt / 6) * (k1_q + 2 * k2_q + 2 * k3_q + k4_q)

    q_new /= np.linalg.norm(q_new)  # Нормализация кватерниона
    return q_new, w_new


def solve_numeric(q0, w0, I_values, r, F, M_ext, t_max, dt):
    t = np.arange(0, t_max, dt)
    q = np.zeros((len(t), 4))
    w = np.zeros((len(t), 3))
    potential_energy = np.zeros(len(t))
    kinetic_energy = np.zeros(len(t))
    total_energy = np.zeros(len(t))

    q[0] = q0
    w[0] = w0

    I = np.diag(I_values)
    I_inv = np.linalg.inv(I)

    for i in range(1, len(t)):
        torque = np.cross(r, F) + M_ext
        q[i], w[i] = rk4_step(q[i - 1], w[i - 1], dt, torque, I, I_inv)

        # Энергии
        rotation_matrix = quaternion_to_rotation_matrix(q[i])
        r_rotated = np.dot(rotation_matrix, r)
        kinetic_energy[i] = compute_kinetic_energy(w[i], I)
        # potential_energy[i] = -np.dot(F, r_rotated)
        potential_energy[i] = -m * g * r_rotated[1]
        total_energy[i] = kinetic_energy[i] + potential_energy[i]

    return t, q, w, potential_energy, kinetic_energy, total_energy


def compute_angular_momentum(w, I):
    return np.dot(I, w)


def compute_kinetic_energy(w, I):
    return 0.5 * np.dot(w, np.dot(I, w))


def update_plot(frame, q, w, ax, I, potential_energy, kinetic_energy, total_energy):
    rotation_matrix = quaternion_to_rotation_matrix(q[frame])

    axes_length = 1
    axes = np.array([
        [axes_length, 0, 0],
        [0, axes_length, 0],
        [0, 0, axes_length],
    ])
    rotated_axes = np.dot(rotation_matrix, axes.T).T

    ax.cla()
    ax.quiver(0, 0, 0, rotated_axes[0, 0], rotated_axes[0, 1], rotated_axes[0, 2], color='r', label='X')
    ax.quiver(0, 0, 0, rotated_axes[1, 0], rotated_axes[1, 1], rotated_axes[1, 2], color='g', label='Y')
    ax.quiver(0, 0, 0, rotated_axes[2, 0], rotated_axes[2, 1], rotated_axes[2, 2], color='b', label='Z')

    ax.set_xlim([-1, 1])
    ax.set_ylim([-1, 1])
    ax.set_zlim([-1, 1])
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_title(f'Frame: {frame}')

    # Расчет кинетического момента
    L = compute_angular_momentum(w[frame], I)
    T = kinetic_energy[frame]
    U = potential_energy[frame]
    E = total_energy[frame]

    # Вывод покоординатного кинетического момента
    print(f"Frame {frame}:")
    print(f"  Angular Momentum (L): X={L[0]:.3f}, Y={L[1]:.3f}, Z={L[2]:.3f}")
    print(f"  Kinetic Energy: {T:.3f}, Potential Energy: {U:.3f}, Total Energy: {E:.3f}")


if __name__ == "__main__":
    q0 = [1, 0, 0, 0]
    w0 = [0, 0, 3]
    I_values = [1, 1, 3]
    r = [0, 0, -1]
    m = 1
    g = -9.8
    F = [0, m * g, 0]
    M_ext = [0, 0, 0]
    t_max = 3
    dt = 0.001

    t, q, w, potential_energy, kinetic_energy, total_energy = solve_numeric(q0, w0, I_values, r, F, M_ext, t_max, dt)

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    ani = FuncAnimation(fig, update_plot, frames=range(0, len(t), 10),
                        fargs=(q, w, ax, np.diag(I_values), potential_energy, kinetic_energy, total_energy),
                        interval=dt * 1000, blit=False)

    ani.save('Lagrang_.gif', writer='imagemagick', fps=30)
    plt.show()
