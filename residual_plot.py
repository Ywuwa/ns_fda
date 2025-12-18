import glob
import matplotlib.pyplot as plt

# Шаблон для поиска файлов
file_pattern = "residual*.txt"

# Список файлов
files = glob.glob(file_pattern)

# Словарь для хранения данных
data = {}

# Чтение данных из файлов
for filename in files:
    with open(filename, "r", encoding="utf-8") as f:
        values = [float(line.strip()) for line in f if line.strip()]
        data[filename] = values

# Построение графика
plt.figure()

for filename, values in data.items():
    x = range(len(values))  # ось времени
    plt.plot(x, values, label=f"k = {filename[-5]}")

# Оформление графика
plt.xlabel("Время")
plt.ylabel("Невязка")
plt.legend()
plt.grid(True)
plt.title("Изменение невязки в зависимости от \n завихренности течения k")

# Сохранение графика в файл
plt.savefig("residual_plot.png", dpi=300, bbox_inches="tight")

# Показать график
plt.show()