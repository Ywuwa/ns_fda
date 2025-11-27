import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401

def read_cube_txt(path):
  """
  Читает файл с "квадратами" разделёнными пустой строкой.
  Возвращает numpy array shape (nz, ny, nx) — порядок Z, Y, X.
  Значения преобразуются в float.
  """
  with open(path, 'r', encoding='utf-8') as f:
      lines = [ln.rstrip() for ln in f]

  blocks = []
  cur_block = []
  for ln in lines:
    if ln.strip() == "":
      if cur_block:
        blocks.append(cur_block)
        cur_block = []
    else:
      # разделяем по пробелам/табуляциям
      row = ln.split()
      cur_block.append(row)
  if cur_block:
    blocks.append(cur_block)

  if not blocks:
    raise ValueError("Файл пуст или некорректен")

  # Поддерживаем различную размерность: 
  # проверяем, что все блоки имеют одинаковую форму
  ny = len(blocks[0])
  nx = len(blocks[0][0])
  nz = len(blocks)
  for b in blocks:
    if len(b) != ny:
      raise ValueError("Несоответствие числа строк в блоках")
    for row in b:
      if len(row) != nx:
        raise ValueError("Несоответствие числа значений в строках")

  arr = np.zeros((nz, ny, nx), dtype=float)
  for k, b in enumerate(blocks):
    for j, row in enumerate(b):
      for i, val in enumerate(row):
        arr[k, j, i] = float(val)
  return arr  # shape (nz, ny, nx)
#==============================================================================

#==============================================================================
def visualize_nodes(values, cube_size=1.0, cmap='viridis', s=5):
  """
  values: numpy array shape (nz, ny, nx)
  cube_size: длина ребра куба (по умолчанию 1.0)
  cmap: matplotlib colormap
  s: размер маркера (шарика)
  """
  nz, ny, nx = values.shape
  s = 500.0 / nx
  # координаты узлов: равномерно от 0 до cube_size вдоль каждой оси
  xs = np.linspace(0, cube_size, nx)
  ys = np.linspace(0, cube_size, ny)
  zs = np.linspace(0, cube_size, nz)

  # создаём сетку координат для всех узлов
  X, Y, Z = np.meshgrid(xs, ys, zs, indexing='xy')  # X,Y shape (ny, nx, nz) -> careful with shapes
  # Перестроим в векторный вид совместимо с values
  # values сейчас (nz, ny, nx) — приведём к (ny, nx, nz) чтобы соответствовать meshgrid indexing='xy'
  vals_reordered = np.transpose(values, (1, 2, 0))
  Xf = X.ravel(); Yf = Y.ravel(); Zf = Z.ravel()
  Vf = vals_reordered.ravel()

  # нормируем цвета по значениям
  norm = plt.Normalize(vmin=np.nanmin(Vf), vmax=np.nanmax(Vf))
  cmap_m = plt.colormaps[cmap]

  fig = plt.figure(figsize=(16, 9))
  ax = fig.add_subplot(111, projection='3d')
  sc = ax.scatter(Xf, Yf, Zf, c=Vf, 
                  alpha=0.5, cmap=cmap_m, norm=norm, s=s, depthshade=True)

  # Настройка осей и видимых границ куба
  ax.set_xlim(0, cube_size); ax.set_ylim(0, cube_size); ax.set_zlim(0, cube_size)
  ax.set_xlabel('X'); ax.set_ylabel('Y'); ax.set_zlabel('Z')
  ax.set_box_aspect((1, 1, 1))  # равные масштабы по осям

  # Добавим цветовую шкалу
  plt.colorbar(sc, ax=ax, pad=0.1, label='Value')

  plt.title('Значения давления в узлах сетки')
  line1 = mlines.Line2D([], [], color='none', label='Custom Text Label')
  line2 = mlines.Line2D([], [], color='none', label='Custom Text Label')
  ax.legend([line1, line2], ('h = 0.1', 't = 0.05 h'), 
            loc='upper left', title='Шаги сетки')
  plt.tight_layout()
  # Получаем изображение из фигуры
  # Сначала рисуем
  fig.canvas.draw()
  # Преобразуем в numpy массив
  frame = fig.canvas.renderer.buffer_rgba()
  #plt.show()
  plt.close(fig) # закрываем фигуру, чтобы не расходовать память
  return frame
#==============================================================================

#========================== PRESSURE TWO SUBPLOTS ============================= 
def visualize_nodes_two_subplots(
    values, values2, cube_size=1.0, cmap='viridis', s=5):
  """
  values: numpy array shape (nz, ny, nx)
  cube_size: длина ребра куба (по умолчанию 1.0)
  cmap: matplotlib colormap
  s: размер маркера (шарика)
  """
  nz, ny, nx = values.shape
  s = 500.0 / nx
  # координаты узлов: равномерно от 0 до cube_size вдоль каждой оси
  xs = np.linspace(0, cube_size, nx)
  ys = np.linspace(0, cube_size, ny)
  zs = np.linspace(0, cube_size, nz)

  # создаём сетку координат для всех узлов
  X, Y, Z = np.meshgrid(xs, ys, zs, indexing='xy')  # X,Y shape (ny, nx, nz) -> careful with shapes
  # Перестроим в векторный вид совместимо с values
  # values сейчас (nz, ny, nx) — приведём к (ny, nx, nz) чтобы соответствовать meshgrid indexing='xy'
  vals_reordered = np.transpose(values, (1, 2, 0))
  Xf = X.ravel(); Yf = Y.ravel(); Zf = Z.ravel()
  Vf = vals_reordered.ravel()

  # нормируем цвета по значениям
  norm = plt.Normalize(vmin=np.nanmin(Vf), vmax=np.nanmax(Vf))
  cmap_m = plt.colormaps[cmap]

  fig = plt.figure(figsize=(16, 9), layout="constrained")
  # first subplot
  ax1 = fig.add_subplot(121, projection='3d')
  sc = ax1.scatter(Xf, Yf, Zf, c=Vf, 
                  alpha=0.5, cmap=cmap_m, norm=norm, s=s, depthshade=True)

  # Настройка осей и видимых границ куба
  ax1.set_xlim(0, cube_size); ax1.set_ylim(0, cube_size); ax1.set_zlim(0, cube_size)
  ax1.set_xlabel('X'); ax1.set_ylabel('Y'); ax1.set_zlabel('Z')
  ax1.set_box_aspect((1, 1, 1))  # равные масштабы по осям
  ax1.set_title('Точное решение')
  
  # second subplot
  vals_reordered = np.transpose(values2, (1, 2, 0))
  Vf = vals_reordered.ravel()
  ax = fig.add_subplot(122, projection='3d')
  sc = ax.scatter(Xf, Yf, Zf, c=Vf, 
                  alpha=0.5, cmap=cmap_m, norm=norm, s=s, depthshade=True)

  # Настройка осей и видимых границ куба
  ax.set_xlim(0, cube_size); ax.set_ylim(0, cube_size); ax.set_zlim(0, cube_size)
  ax.set_xlabel('X'); ax.set_ylabel('Y'); ax.set_zlabel('Z')
  ax.set_box_aspect((1, 1, 1))  # равные масштабы по осям
  ax.set_title('Решение системы')

  # colorbar, legend, title
  fig.suptitle('Значения давления в узлах сетки', fontsize=16)
  # Добавим цветовую шкалу
  #plt.colorbar(sc, ax=ax, pad=0.1, label='Value')
  sm = plt.cm.ScalarMappable(cmap=cmap_m, norm=norm)
  sm.set_array([])
  fig.colorbar(sm, ax=[ax1, ax], pad=0.1, label='Модуль скорости')

  line1 = mlines.Line2D([], [], color='none', label='Custom Text Label')
  line2 = mlines.Line2D([], [], color='none', label='Custom Text Label')
  fig.legend([line1, line2], ('h = 0.1', 't = 0.05 h'), 
            loc='upper left', title='Шаги сетки')

  # Получаем изображение из фигуры
  # Сначала рисуем
  fig.canvas.draw()
  # Преобразуем в numpy массив
  frame = fig.canvas.renderer.buffer_rgba()
  #plt.show()
  plt.close(fig) # закрываем фигуру, чтобы не расходовать память
  return frame
#==============================================================================

#==============================================================================
def visualize_arrows(v1, v2, v3,
                     cube_size=1.0, step_index=10, scale=0.2, cmap='plasma'):
  """
  values: shape (nz, ny, nx)
  step_index: берём узлы с индексами кратными этому числу.
  scale: визуальный множитель для длины стрелок.
  """
  nz, ny, nx = v1.shape
  scale = 2.0 / nx
  step_index = nx // 10
  xs = np.linspace(0, cube_size, nx)
  ys = np.linspace(0, cube_size, ny)
  zs = np.linspace(0, cube_size, nz)

  idx_i = list(range(0, nx, step_index))
  idx_j = list(range(0, ny, step_index))
  idx_k = list(range(0, nz, step_index))

  Xf = []; Yf = []; Zf = []
  U = []; V = []; W = []

  for k in idx_k:
    for j in idx_j:
      for i in idx_i:
        x = xs[i]; y = ys[j]; z = zs[k]
        Xf.append(x); Yf.append(y); Zf.append(z)
        ux = v1[k, j, i]; uy = v2[k, j, i]; uz = v3[k, j, i]
        U.append(ux); V.append(uy); W.append(uz)

  if len(Xf) == 0:
    raise ValueError("Нет точек для отображения. Увеличьте размер сетки или уменьшите step_index.")

  U = np.array(U); V = np.array(V); W = np.array(W)
  magn = np.sqrt(U**2 + V**2 + W**2)
  maxm = magn.max() if magn.max() != 0 else 1.0
  Un = U / maxm; Vn = V / maxm; Wn = W / maxm

  U_plot = Un * scale * cube_size
  V_plot = Vn * scale * cube_size
  W_plot = Wn * scale * cube_size

  C = np.array(magn)
  norm = plt.Normalize(vmin=np.nanmin(C), vmax=np.nanmax(C))
  cmap_m = plt.colormaps[cmap]

  fig = plt.figure(figsize=(16,9))
  ax = fig.add_subplot(111, projection='3d')

  ax.quiver(np.array(Xf), np.array(Yf), np.array(Zf),
            U_plot, V_plot, W_plot,
            length=1.0, normalize=False, colors=cmap_m(norm(C)), linewidths=0.7)

  ax.set_xlim(0, cube_size); ax.set_ylim(0, cube_size); ax.set_zlim(0, cube_size)
  ax.set_xlabel('X'); ax.set_ylabel('Y'); ax.set_zlabel('Z')
  ax.set_box_aspect((1,1,1))

  sm = plt.cm.ScalarMappable(cmap=cmap_m, norm=norm)
  sm.set_array([])
  plt.colorbar(sm, ax=ax, pad=0.1, label='Модуль скорости')

  plt.title('Поле скоростей')
  line1 = mlines.Line2D([], [], color='none', label='Custom Text Label')
  line2 = mlines.Line2D([], [], color='none', label='Custom Text Label')
  ax.legend([line1, line2], ('h = 0.1', 't = 0.05 h'), 
            loc='upper left', title='Шаги сетки')

  plt.tight_layout()
  #plt.show()
  # Получаем изображение из фигуры
  # Сначала рисуем
  fig.canvas.draw()
  # Преобразуем в numpy массив
  frame = fig.canvas.renderer.buffer_rgba()
  #plt.show()
  plt.close(fig) # закрываем фигуру, чтобы не расходовать память
  return frame
#==============================================================================

#==============================================================================
def visualize_arrows_two_subplots(v1, v2, v3, v1_2, v2_2, v3_2,
                     cube_size=1.0, step_index=10, scale=0.2, cmap='plasma'):
  """
  values: shape (nz, ny, nx)
  step_index: берём узлы с индексами кратными этому числу.
  scale: визуальный множитель для длины стрелок.
  """
  nz, ny, nx = v1.shape
  scale = 2.0 / nx
  step_index = nx // 10
  xs = np.linspace(0, cube_size, nx)
  ys = np.linspace(0, cube_size, ny)
  zs = np.linspace(0, cube_size, nz)

  idx_i = list(range(0, nx, step_index))
  idx_j = list(range(0, ny, step_index))
  idx_k = list(range(0, nz, step_index))

  Xf = []; Yf = []; Zf = []
  U = []; V = []; W = []
  U_2 = []; V_2 = []; W_2 = []

  for k in idx_k:
    for j in idx_j:
      for i in idx_i:
        x = xs[i]; y = ys[j]; z = zs[k]
        Xf.append(x); Yf.append(y); Zf.append(z)
        
        ux = v1[k, j, i]; uy = v2[k, j, i]; uz = v3[k, j, i]
        U.append(ux); V.append(uy); W.append(uz)
        
        ux = v1_2[k, j, i]; uy = v2_2[k, j, i]; uz = v3_2[k, j, i]
        U_2.append(ux); V_2.append(uy); W_2.append(uz)

  if len(Xf) == 0:
    raise ValueError("Нет точек для отображения. Увеличьте размер сетки или уменьшите step_index.")

  U = np.array(U); V = np.array(V); W = np.array(W)
  magn = np.sqrt(U**2 + V**2 + W**2)
  maxm = magn.max() if magn.max() != 0 else 1.0
  Un = U / maxm; Vn = V / maxm; Wn = W / maxm

  U_plot = Un * scale * cube_size
  V_plot = Vn * scale * cube_size
  W_plot = Wn * scale * cube_size

  C = np.array(magn)
  norm = plt.Normalize(vmin=np.nanmin(C), vmax=np.nanmax(C))
  cmap_m = plt.colormaps[cmap]

  fig = plt.figure(figsize=(16,9), layout="constrained")
  # first plot
  ax1 = fig.add_subplot(121, projection='3d')

  ax1.quiver(np.array(Xf), np.array(Yf), np.array(Zf),
            U_plot, V_plot, W_plot,
            length=1.0, normalize=False, colors=cmap_m(norm(C)), linewidths=0.7)

  ax1.set_xlim(0, cube_size); ax1.set_ylim(0, cube_size); ax1.set_zlim(0, cube_size)
  ax1.set_xlabel('X'); ax1.set_ylabel('Y'); ax1.set_zlabel('Z')
  ax1.set_box_aspect((1,1,1))
  ax1.set_title('Точное решение')
  
  # second plot
  #magn = np.sqrt(U_2**2 + V_2**2 + W_2**2)
  magn = np.sqrt(U**2 + V**2 + W**2)
  maxm = magn.max() if magn.max() != 0 else 1.0
  Un = U_2 / maxm; Vn = V_2 / maxm; Wn = W_2 / maxm

  U_plot = Un * scale * cube_size
  V_plot = Vn * scale * cube_size
  W_plot = Wn * scale * cube_size
  C = np.array(magn)
  norm = plt.Normalize(vmin=np.nanmin(C), vmax=np.nanmax(C))
  cmap_m = plt.colormaps[cmap]
  ax = fig.add_subplot(122, projection='3d')

  ax.quiver(np.array(Xf), np.array(Yf), np.array(Zf),
            U_plot, V_plot, W_plot,
            length=1.0, normalize=False, colors=cmap_m(norm(C)), linewidths=0.7)

  ax.set_xlim(0, cube_size); ax.set_ylim(0, cube_size); ax.set_zlim(0, cube_size)
  ax.set_xlabel('X'); ax.set_ylabel('Y'); ax.set_zlabel('Z')
  ax.set_box_aspect((1,1,1))
  ax.set_title('Решение системы')
  
  # colorbar, legend, title
  fig.suptitle('Поле скоростей', fontsize=16)
  sm = plt.cm.ScalarMappable(cmap=cmap_m, norm=norm)
  sm.set_array([])
  fig.colorbar(sm, ax=[ax1, ax], pad=0.1, label='Модуль скорости')
  #fig.subplots_adjust(right=0.8)
  #cbar_ax = fig.add_axes([0.9, 0.2, 0.02, 0.6])
  #fig.colorbar(sm, cax=cbar_ax)

  line1 = mlines.Line2D([], [], color='none', label='Custom Text Label')
  line2 = mlines.Line2D([], [], color='none', label='Custom Text Label')
  fig.legend([line1, line2], ('h = 0.1', 't = 0.05 h'), 
            loc='upper left', title='Шаги сетки')

  #plt.tight_layout(pad=0)
  # Получаем изображение из фигуры
  # Сначала рисуем
  fig.canvas.draw()
  # Преобразуем в numpy массив
  frame = fig.canvas.renderer.buffer_rgba()
  #plt.show()
  plt.close(fig) # закрываем фигуру, чтобы не расходовать память
  return frame
#==============================================================================

#=================================== MAIN =====================================
if __name__ == '__main__':
  # Пример использования:
  # filename = 'data.txt'
  # filename = input("Путь к файлу с данными: ").strip()
  num = 10
  filename = "helical_p/" + str(num) + ".txt"
  filename2 = "p/" + str(num) + ".txt"
  values = read_cube_txt(filename)  # shape (nz, ny, nx)
  values2 = read_cube_txt(filename2)  # shape (nz, ny, nx)
  visualize_nodes(values, cube_size=1.0, cmap='viridis', s=80)
  visualize_nodes_two_subplots(values, values2)
  
  v1_file = "helical_v1/" + str(num) + ".txt"
  v2_file = "helical_v2/" + str(num) + ".txt"
  v3_file = "helical_v3/" + str(num) + ".txt"
  v1_val = read_cube_txt(v1_file)
  v2_val = read_cube_txt(v2_file)
  v3_val = read_cube_txt(v3_file)
  
  v1 = "v1/" + str(num) + ".txt";
  v2 = "v2/" + str(num) + ".txt"
  v3 = "v3/" + str(num) + ".txt"
  v1_val_2 = read_cube_txt(v1)
  v2_val_2 = read_cube_txt(v2) 
  v3_val_2 = read_cube_txt(v3)
  visualize_arrows_two_subplots(v1_val, v2_val, v3_val, v1_val_2, v2_val_2, v3_val_2)