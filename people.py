from config import INFO
import subprocess
import time

def wait_a_minute():
  time.sleep(60)

def is_working():
  return bool(subprocess.Popen('top -u zhangz -b -n 1 | grep python', shell=True, stdout=subprocess.PIPE).communicate()[0])

def wait_until_idle():
  if is_working():
    INFO('waiting indefinitely for other jobs to finish...')
  while is_working():
    wait_a_minute()

# vim: et sw=2
