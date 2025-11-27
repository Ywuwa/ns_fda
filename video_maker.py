import matplotlib.pyplot as plt
import matplotlib.animation as animation
import multiprocessing

# import module
#==============================================================================
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
#======================== FRAME WORK FOR PRESSURE =============================
def get_frame_pressure_pair(pfile, pfile2):
  values = visualize.read_cube_txt(pfile)  # shape (nz, ny, nx)
  values2 = visualize.read_cube_txt(pfile2)  # shape (nz, ny, nx)
  frame = visualize.visualize_nodes_two_subplots(values, values2)
  return frame
#==============================================================================

#======================== FRAME WORK FOR VELOCITY =============================
def get_frame_velocity_pair(v1_file, v2_file, v3_file,
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

def create_video_velocity_pair(q, file, video):
  num_frames = 300
  fps = 5
  fig, ax = plt.subplots()
  ax.set_axis_off()  # отключим оси, чтобы было похоже на видео
  
  # Создаем объект для записи видео через matplotlib
  Writer = animation.FFMpegWriter
  writer = Writer(fps=fps, metadata=dict(artist='Me'), bitrate=1800)
  
  with writer.saving(fig, video, dpi=300):
    for i in range(num_frames):
      v1_file = "helical_v1/" + f"{i}" + ".txt"
      v2_file = "helical_v2/" + f"{i}" + ".txt"
      v3_file = "helical_v3/" + f"{i}" + ".txt"
      
      v1_f2 = "v1/" + f"{i}" + ".txt"
      v2_f2 = "v2/" + f"{i}" + ".txt"
      v3_f2 = "v3/" + f"{i}" + ".txt"
      frame = get_frame_velocity_pair(
        v1_file, v2_file, v3_file, v1_f2, v2_f2, v3_f2)
      # Отобразим кадр в `imshow`
      ax.imshow(frame)
      # Добавляем следующий кадр
      writer.grab_frame()
      # Очистим фигуру для следующего кадра
      ax.clear()
      ax.set_axis_off()
      if (i % 10 == 0):
        q.put(f"{i} of {num_frames} for {video}")
  
  q.put(f"Видео сохранено как {video}.mp4")
#==============================================================================

#============================ CREATE VIDEO ====================================
def create_video_pressure_pair(q, file, video):
  num_frames = 300
  fps = 12
  fig, ax = plt.subplots()
  ax.set_axis_off()  # отключим оси, чтобы было похоже на видео
  
  # Создаем объект для записи видео через matplotlib
  Writer = animation.FFMpegWriter
  writer = Writer(fps=fps, metadata=dict(artist='Me'), bitrate=1800)

  with writer.saving(fig, video, dpi=300):
    for i in range(num_frames):
      pfile = "helical_p/" + f"{i}" + ".txt"
      pfile2 = "p/" + f"{i}" + ".txt"
      frame = get_frame_pressure_pair(pfile, pfile2)
      # Отобразим кадр в `imshow`
      ax.imshow(frame)
      # Добавляем следующий кадр
      writer.grab_frame()
      # Очистим фигуру для следующего кадра
      ax.clear()
      ax.set_axis_off()
      if (i % 10 == 0):
        q.put(f"{i} of {num_frames} for {video}")
  
  q.put(f"Видео сохранено как {video}.mp4")
#=================================== MAIN =====================================
if __name__ == "__main__":
  queue = multiprocessing.Queue()
  # Определяем файлы и имена выводимых видео
  file1 = 'python_logger.txt'  # заменить на путь к первому файлу
  file2 = 'python_logger.txt'  # заменить на путь ко второму файлу

  video1 = 'pressure_video.mp4'
  video2 = 'velocity_video.mp4'

  # Создаём процессы для каждого видео
  p1 = multiprocessing.Process(
    target=create_video_pressure_pair, args=(queue, file1, video1))
  p2 = multiprocessing.Process(
    target=create_video_velocity_pair, args=(queue, file2, video2))

  # Запускаем процессы
  p1.start()
  p2.start()

  # Ожидаем завершения
  p1.join()
  p2.join()
  while not queue.empty():
    print(queue.get())
  print("Оба видео созданы.")
#===================================================================================