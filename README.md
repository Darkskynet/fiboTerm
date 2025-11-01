# 🌌 fiboTerm

**fiboTerm** is a live, terminal-based **Fibonacci / Golden Spiral visualizer** — a hypnotic piece of mathematical art written entirely in Python.

It draws a single continuous golden spiral that grows outward from the center of your terminal using smooth ASCII or high-resolution **Unicode Braille** characters.  
Colors, pulsing, trails, and even spinning motion are fully configurable.

---

## ✨ Features

- 🌀 **Single-arm Golden Spiral** (logarithmic spiral where each 90° grows by the golden ratio φ)
- 🟣 **High-Resolution Braille Mode** (`--braille`) for ultra-smooth curves  
- 🌈 **Full ANSI Color Support** (`--color rainbow | pulse | none`)
- 💫 **Trail, Spin, and Pulse** animation modes  
- 🔁 **Hypnotic Mode** (`--hypnotic`) for infinite looping
- 🧮 Optional on-screen **Fibonacci sequence** and **φ (phi)** display
- 🧠 **Auto-scales** to your terminal window size

No external libraries — just Python’s standard library.

---

## 🚀 Usage

```bash
git clone https://github.com/Darkskynet/fiboTerm.git
cd fiboTerm
python fiboTerm.py [options]
````

### Examples

#### Classic rainbow ASCII spiral

```bash
python fiboTerm.py --color rainbow
```

#### High-resolution Braille mode

```bash
python fiboTerm.py --braille --color rainbow --turns 10
```

#### Hypnotic spinning mode

```bash
python fiboTerm.py --braille --trail --spin --hypnotic
```

#### Minimalist monochrome

```bash
python fiboTerm.py --color none
```

---

## ⚙️ Options

| Flag             | Description                                        |
| ---------------- | -------------------------------------------------- |
| `--braille`      | Use Unicode Braille grid (2×4 subpixels per cell). |
| `--trail`        | Enable fading trails behind the spiral.            |
| `--trail-decay`  | Decay rate of the trail (default `0.02`).          |
| `--spin`         | Rotate the completed spiral in place.              |
| `--spin-steps`   | Number of frames per spin cycle (default `80`).    |
| `--spin-delay`   | Delay between spin frames (default `0.02s`).       |
| `--pulse`        | Add rhythmic brightness pulsing.                   |
| `--color`        | Color mode: `rainbow`, `pulse`, or `none`.         |
| `--color-source` | Color change source: `time` or `distance`.         |
| `--goldenratio`  | Show the golden ratio φ ≈ 1.61803.                 |
| `--fibs`         | Show Fibonacci numbers as they appear.             |
| `--hypnotic`     | Loop forever for an infinite spiral.               |
| `--turns`        | Number of 90° turns to draw (default `6`).         |
| `--delay`        | Delay between points in seconds.                   |
| `--density`      | Points per degree of curvature (default `4`).      |
| `--batch`        | Render every N points to balance smoothness/speed. |

---

## 📚 Mathematical Background

The **Fibonacci sequence** is a series of numbers where each term is the sum of the previous two:

> 0, 1, 1, 2, 3, 5, 8, 13, 21, 34, ...

As `n` grows, the ratio between successive Fibonacci numbers approaches the **golden ratio**:

> φ = (1 + √5) / 2 ≈ 1.61803

The **golden spiral** is a special type of **logarithmic spiral** defined by:

> r(θ) = a × φ^(θ / 90°)

Every quarter turn (90°) increases the radius by a factor of φ.

### Further Reading

* [Wikipedia: Fibonacci number](https://en.wikipedia.org/wiki/Fibonacci_number)
* [Wikipedia: Golden ratio](https://en.wikipedia.org/wiki/Golden_ratio)
* [Wikipedia: Golden spiral](https://en.wikipedia.org/wiki/Golden_spiral)
* [Wikipedia: Logarithmic spiral](https://en.wikipedia.org/wiki/Logarithmic_spiral)
* [Unicode Braille Patterns (U+2800–U+28FF)](https://en.wikipedia.org/wiki/Braille_Patterns)

---

## 🧩 Technical Notes

* Written for **Python 3.8+**
* Works best in a 256-color capable terminal (Windows Terminal, iTerm2, Alacritty, Kitty, etc.)
* No dependencies outside the standard library
* Automatically adjusts for your terminal’s width/height
* Cleans up cursor visibility and color state on exit

---

## 🖼️ Screenshots

<img width="814.5" height="610" alt="image" src="https://github.com/user-attachments/assets/9fafce18-92bb-4ff9-9370-9b191ba39e4e" />
<img width="408" height="596" alt="image" src="https://github.com/user-attachments/assets/a05d4f6d-009f-4ec8-adcc-479442bb3255" />



---

## 💡 Ideas for Future Enhancements

* Keyboard-interactive controls (change color mode live)
* Sound synthesis linked to radius or color
* PNG export via Pillow
* Text backgrounds or ambient ASCII art
* Smooth fade-in/out transitions between modes

---

## 🧑‍💻 Author

Created by **Darkskynet**

---

## 🪐 License

MIT License – free to use, modify, and remix.

---

> “Mathematics is the poetry of logical ideas.” – Albert Einstein

---
