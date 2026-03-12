import matplotlib.pyplot as plt
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

        x, y, _ = map(float, re.findall(r"[-+]?\d*\.?\d+", lines[i]))
        i += 1

        largura = float(re.findall(r"\d+", lines[i])[0])
        i += 1

        comprimento = float(re.findall(r"\d+", lines[i])[0])
        i += 1

        i += 1  # Altura (ignored)

        items.append({
            "id": item_id,
            "x": x,
            "y": y,
            "w": largura,
            "h": comprimento
        })

    return items


def plot_xy(items):
    fig, ax = plt.subplots(figsize=(8, 6))

    for it in items:
        rect = plt.Rectangle(
            (it["x"], it["y"]),
            it["w"],
            it["h"],
            fill=False,
            linewidth=2
        )
        ax.add_patch(rect)

        ax.text(
            it["x"] + it["w"] / 2,
            it["y"] + it["h"] / 2,
            f'ID {it["id"]}',
            ha="center",
            va="center"
        )

    # 🔥 THIS IS THE KEY PART
    ax.relim()
    ax.autoscale_view()

    ax.set_xlabel("X (Largura)")
    ax.set_ylabel("Y (Comprimento)")
    ax.set_aspect("equal", adjustable="box")
    ax.grid(True)

    plt.show()


if __name__ == "__main__":
    items = read_items("items.txt")
    print(items)  # debug: should NOT be empty
    plot_xy(items)

