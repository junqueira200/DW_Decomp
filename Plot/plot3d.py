import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import re


def read_items(filename):
    with open(filename, "r") as f:
        lines = [line.strip() for line in f if line.strip()]

    n = int(lines[0])
    items = []

    i = 1
    for _ in range(n):
        item_id = int(re.search(r"\d+", lines[i]).group())
        i += 1

        x, y, z = map(float, re.findall(r"[-+]?\d*\.?\d+", lines[i]))
        i += 1

        largura = float(re.findall(r"\d+", lines[i])[0])
        i += 1

        comprimento = float(re.findall(r"\d+", lines[i])[0])
        i += 1

        altura = float(re.findall(r"\d+", lines[i])[0])
        i += 1

        items.append({
            "id": item_id,
            "x": x,
            "y": y,
            "z": z,
            "dx": largura,
            "dy": comprimento,
            "dz": altura
        })

    return items


def draw_box(ax, x, y, z, dx, dy, dz):
    vertices = [
        [(x, y, z), (x+dx, y, z), (x+dx, y+dy, z), (x, y+dy, z)],
        [(x, y, z+dz), (x+dx, y, z+dz), (x+dx, y+dy, z+dz), (x, y+dy, z+dz)],
        [(x, y, z), (x+dx, y, z), (x+dx, y, z+dz), (x, y, z+dz)],
        [(x, y+dy, z), (x+dx, y+dy, z), (x+dx, y+dy, z+dz), (x, y+dy, z+dz)],
        [(x, y, z), (x, y+dy, z), (x, y+dy, z+dz), (x, y, z+dz)],
        [(x+dx, y, z), (x+dx, y+dy, z), (x+dx, y+dy, z+dz), (x+dx, y, z+dz)]
    ]

    ax.add_collection3d(Poly3DCollection(vertices, alpha=0.4))


def plot_items(items):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection="3d")

    for it in items:
        draw_box(ax, it["x"], it["y"], it["z"],
                 it["dx"], it["dy"], it["dz"])

        ax.text(it["x"], it["y"], it["z"],
                f'ID {it["id"]}', fontsize=9)

    ax.set_xlabel("X (Largura)")
    ax.set_ylabel("Y (Comprimento)")
    ax.set_zlabel("Z (Altura)")

    ax.set_box_aspect([1, 1, 1])
    plt.show()


if __name__ == "__main__":
    items = read_items("items.txt")
    plot_items(items)

