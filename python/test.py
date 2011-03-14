import Image
import ImageDraw

x = Image.new('RGB', (500,500))
draw = ImageDraw.Draw(x)
draw.ellipse((0,0,100,100), fill='red')
x.save('out.png', 'png')
