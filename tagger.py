class TagCollection:
  def __init__(self, taglist):
    pass

class Tagger:
  def __init__(self, image, info):
    self.image = image
    self.info = info

  def get_frustrum(self):
    pass # TODO

  def get_distance(self, tag):
    pass # TODO

  def get_pixels(self):
    pass # TODO

  def get_tag_points(self):
    "Returns collection of (tag, pixel) pairs"
    DIST_THRESHOLD = 1.0 # meters
    possible_tags = self.get_frustrum()
    tag_positions = {}
    for pixel in self.get_pixels():
      for tag in possible_tags:
        dist = self.get_distance(tag)
        if dist < DIST_THRESHOLD:
          tag_positions[tag] = (dist, pixel)
    for tag in tag_positions:
      places = tag_positions[tag]
      if places:
        yield (tag, min(places)[1])

  def tag(self, output):
    points = self.get_tag_points()
    # TODO write points to output

# vim: et sw=2
