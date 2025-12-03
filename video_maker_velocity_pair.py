import matplotlib.pyplot as plt
import matplotlib.animation as animation

try:
    # Trying to find module in the parent package
    from . import visualize
    print(True)
    
except ImportError:
    print('Relative import failed')

try:
    # Trying to find module on sys.path
    import visualize
    print(True)
except ModuleNotFoundError:
    print('Absolute import failed')

def get_frame(v1_file, v2_file, v3_file,
              v1_f2, v2_f2, v3_f2):
  v1_val = visualize.read_cube_txt(v1_file)
  v2_val = visualize.read_cube_txt(v2_file)
  v3_val = visualize.read_cube_txt(v3_file) # shape (nz, ny, nx)
  v1_val2 = visualize.read_cube_txt(v1_f2)
  v2_val2 = visualize.read_cube_txt(v2_f2)
  v3_val2 = visualize.read_cube_txt(v3_f2) # shape (nz, ny, nx)
  frame = visualize.visualize_arrows_two_subplots(
    v1_val, v2_val, v3_val, v1_val2, v2_val2, v3_val2)
  return frame

num_frames = 330
fps = 5
fig, ax = plt.subplots()
ax.set_axis_off()  # отключим оси, чтобы было похоже на видео

# Создаем объект для записи видео через matplotlib
Writer = animation.FFMpegWriter
writer = Writer(fps=fps, metadata=dict(artist='Me'), bitrate=1800)

with writer.saving(fig, "velocity_video.mp4", dpi=300):
  for i in range(num_frames):
    v1_file = "helical_v1/" + f"{i}" + ".txt"
    v2_file = "helical_v2/" + f"{i}" + ".txt"
    v3_file = "helical_v3/" + f"{i}" + ".txt"
    
    v1_f2 = "v1/" + f"{i}" + ".txt"
    v2_f2 = "v2/" + f"{i}" + ".txt"
    v3_f2 = "v3/" + f"{i}" + ".txt"
    frame = get_frame(v1_file, v2_file, v3_file, v1_f2, v2_f2, v3_f2)
    # Отобразим кадр в `imshow`
    ax.imshow(frame)
    # Добавляем следующий кадр
    writer.grab_frame()
    # Очистим фигуру для следующего кадра
    ax.clear()
    ax.set_axis_off()
    if (i % 5 == 0):
      print(f"{i} of {num_frames}")

print("Видео сохранено как velocity_video.mp4")