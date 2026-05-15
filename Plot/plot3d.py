import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import re


# Truck dimensions
TRUCK_WIDTH = 13500
TRUCK_LENGTH = 2450
TRUCK_HEIGHT = 2700


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
            "dx": largura,      # width
            "dy": comprimento,  # length
            "dz": altura        # height
        })

    return items


def draw_box(ax, x, y, z, dx, dy, dz,
             alpha=0.4, facecolor=None, edgecolor="black"):

    vertices = [
        [(x, y, z), (x+dx, y, z), (x+dx, y+dy, z), (x, y+dy, z)],
        [(x, y, z+dz), (x+dx, y, z+dz), (x+dx, y+dy, z+dz), (x, y+dy, z+dz)],
        [(x, y, z), (x+dx, y, z), (x+dx, y, z+dz), (x, y, z+dz)],
        [(x, y+dy, z), (x+dx, y+dy, z), (x+dx, y+dy, z+dz), (x, y+dy, z+dz)],
        [(x, y, z), (x, y+dy, z), (x, y+dy, z+dz), (x, y, z+dz)],
        [(x+dx, y, z), (x+dx, y+dy, z), (x+dx, y+dy, z+dz), (x+dx, y, z+dz)]
    ]

    poly = Poly3DCollection(
        vertices,
        alpha=alpha,
        edgecolors=edgecolor,
        linewidths=0.5,
        facecolors=facecolor
    )

    ax.add_collection3d(poly)


def plot_items(items):
    fig = plt.figure(figsize=(14, 8))
    ax = fig.add_subplot(111, projection="3d")

    # Truck boundary
    draw_box(
        ax,
        0, 0, 0,
        TRUCK_WIDTH,
        TRUCK_LENGTH,
        TRUCK_HEIGHT,
        alpha=0.02,
        edgecolor="red"
    )

    for it in items:

        facecolor=None
        if it["id"] == 65:
            facecolor='cyan'

        draw_box(
            ax,
            it["x"],
            it["y"],
            it["z"],
            it["dx"],   # largura
            it["dy"],   # comprimento
            it["dz"],    # altura
            facecolor=facecolor
        )

        ax.text(
            it["x"],
            it["y"],
            it["z"],
            f'ID {it["id"]}',
            fontsize=9
        )

    ax.set_xlabel("X (Largura)")
    ax.set_ylabel("Y (Comprimento)")
    ax.set_zlabel("Z (Altura)")

    ax.set_xlim(0, TRUCK_WIDTH)
    ax.set_ylim(0, TRUCK_LENGTH)
    ax.set_zlim(0, TRUCK_HEIGHT)

    ax.set_box_aspect([
        TRUCK_WIDTH,
        TRUCK_LENGTH,
        TRUCK_HEIGHT
    ])

    plt.show()


if __name__ == "__main__":
    items = read_items("items.txt")
    plot_items(items)
