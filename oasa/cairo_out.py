#--------------------------------------------------------------------------
#     This file is part of OASA - a free chemical python library
#     Copyright (C) 2003-2008 Beda Kosata <beda@zirael.org>

#     This program is free software; you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#     the Free Software Foundation; either version 2 of the License, or
#     (at your option) any later version.

#     This program is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.

#     Complete text of GNU GPL can be found in the file gpl.txt in the
#     main directory of the program

#--------------------------------------------------------------------------


import cairo
import geometry
import math
import misc
import operator
import copy
import sys


class cairo_out:

  """
  This object is used to draw OASA molecules using cairo. Cairo supports different
  'surfaces' which represent different file formats.
  This object implements PNG file drawing, but should be general enough to work with other
  formats, provided modified version of create_surface and write_surface methods are provided
  when this class is subclassed.

  Usage:

  # attributes can be set in constructor or afterwards as normal attributes
  # default options are set ans commented in the default_options dictionary below
  c = cairo_out( scaling=2.0, margin=10, font_size=20, bond_width=6)
  c.show_hydrogens_on_hetero = True
  c.mol_to_cairo( mol, 'outfile.png') # mol is oasa molecule
  
  """

  _temp_margin = 200

  # all the following values are settable using the constructor, e.g.
  # cairo_out( margin=20, bond_width=3.0)
  # all metrics is scaled properly, the values corespond to pixels only
  # when scaling is 1.0
  default_options = {
    'scaling': 1.0,
    'show_hydrogens_on_hetero': False,
    'margin': 15,
    'line_width': 2.0,
    # how far second bond is drawn
    'bond_width': 6.0,
    'wedge_width': 6.0,
    'font_name': "Arial",
    'font_size': 20,
    'background_color': (1,1,1),
    'color_atoms': True,
    'color_bonds': True,
    'space_around_atom': 2,
    # how much to shorten second line of double and triple bonds (0-0.5)
    'bond_second_line_shortening': 0.15,
    # the following two are just for playing
    # - without antialiasing the output is ugly
    'antialias_text': True,
    'antialias_drawing': True,
    # this will only change appearance of overlapping text
    'add_background_to_text': False,
    }

  atom_colors = {'O': (1,0,0),
                 'N': (0,0,1),
                 'S': (.5,.5,0),
                 'Cl': (0,0.8,0),
                 'Br': (.5,0,0),
                 }


  def __init__( self, **kw):
    for k,v in self.__class__.default_options.iteritems():
      setattr( self, k, v)
    # list of paths that contribute to the bounding box (probably no edges)
    self._vertex_to_bbox = {} # vertex-to-bbox mapping
    self._bboxes = [] # for overall bbox calcualtion
    for k,v in kw.items():
      if k in self.__class__.default_options:
        setattr( self, k, v)
      else:
        raise Exception( "unknown attribute '%s' passed to constructor" % k)


  def draw_mol( self, mol):
    if not self.surface:
      raise Exception( "You must initialize cairo surface before drawing, use 'create_surface' to do it.")
    self.molecule = mol
    for v in mol.vertices:
      self._draw_vertex( v)
    for e in copy.copy( mol.edges):
      self._draw_edge( e)

  def create_surface( self, w, h):
    """currently implements PNG writting, but might be overriden to write other types;
    w and h are minimal estimated width and height"""
    # trick - we use bigger surface and then copy from it to a new surface and crop
    self.surface = cairo.ImageSurface( cairo.FORMAT_ARGB32, w, h)


  def write_surface( self):
    """currently implements PNG writting, but might be overriden to write other types"""
    # Because it is not possible to calculate the bounding box of a drawing before its drawn (mainly
    # because we don't know the size of text items), this object internally draws to a surface with
    # large margins and the saves a cropped version into a file (this is done by drawing to a new
    # surface and using the old one as source.
    #self.surface.write_to_png( self.filename)
    #return 
    # real width and height
    x1, y1, x2, y2 = self._get_bbox()
    x1, y1 = self.context.user_to_device( x1, y1)
    x2, y2 = self.context.user_to_device( x2, y2)
    width = int( x2 - x1 + 2*self.margin*self.scaling)
    height = int( y2 - y1 + 2*self.margin*self.scaling)
    surface = cairo.ImageSurface( cairo.FORMAT_ARGB32, width, height)
    context = cairo.Context( surface)
    context.set_source_rgb( *self.background_color)
    context.rectangle( 0, 0, width, height)
    context.fill()
    context.set_source_surface( self.surface, -x1+self.margin*self.scaling, -y1+self.margin*self.scaling)
    context.rectangle( 0, 0, width, height)
    context.fill()
    context.show_page()
    surface.write_to_png( self.filename)
    surface.finish()


  def mol_to_cairo( self, mol, filename):
    x1, y1, x2, y2 = None, None, None, None
    for v in mol.vertices:
      v.y = -v.y # flip coords - molfiles have them the other way around
      if x1 == None or x1 > v.x:
        x1 = v.x
      if x2 == None or x2 < v.x:
        x2 = v.x
      if y1 == None or y1 > v.y:
        y1 = v.y
      if y2 == None or y2 < v.y:
        y2 = v.y
    w = int( x2 - x1)
    h = int( y2 - y1)
    self._bboxes.append( (x1,y1,x2,y2))

    self.filename = filename

    # create cairo surface for drawing
    _w = int( w+2*self.scaling*self._temp_margin)
    _h = int( h+2*self.scaling*self._temp_margin)
    self.create_surface( _w, _h)

    self.context = cairo.Context( self.surface)
    if not self.antialias_drawing:
      self.context.set_antialias( cairo.ANTIALIAS_NONE)
    if not self.antialias_text:
      options = self.context.get_font_options()
      options.set_antialias( cairo.ANTIALIAS_NONE)
      self.context.set_font_options( options)
    self.context.translate( -x1*self.scaling+self.scaling*self._temp_margin, -y1*self.scaling+self.scaling*self._temp_margin)
    self.context.scale( self.scaling, self.scaling)
    self.context.rectangle( x1, y1, w, h)
    self.context.new_path()
    self.context.set_source_rgb( 0, 0, 0)
    self.draw_mol( mol)
    self.context.show_page()
    # write the content to the file
    self.write_surface()
    self.surface.finish()
    # flip y coordinates back
    for v in mol.vertices:
      v.y = -v.y


  def _draw_edge( self, e):
    def draw_plain_or_colored_line( _start, _end, second=False):
      """second means if this is not the main line, drawing might be different"""
      if not has_shown_vertex:
        if not second:
          self._draw_line( _start, _end, line_width=self.line_width, capstyle=cairo.LINE_CAP_ROUND)
        else:
          self._draw_line( _start, _end, line_width=self.line_width, capstyle=cairo.LINE_CAP_BUTT)
      else:
        self._draw_colored_line( _start, _end, line_width=self.line_width, start_color=color1, end_color=color2)

    def draw_plain_or_colored_wedge( _start, _end):
      x1, y1 = _start
      x2, y2 = _end
      x, y, x0, y0 = geometry.find_parallel( x1, y1, x2, y2, self.wedge_width/2.0)
      xa, ya, xb, yb = geometry.find_parallel( x1, y1, x2, y2, self.line_width/2.0) 
      # no coloring now
      if not has_shown_vertex:
        self._create_cairo_path( [(xa, ya), (x0, y0), (2*x2-x0, 2*y2-y0), (2*x1-xa, 2*y1-ya)], closed=True)
        self.context.set_source_rgb( 0,0,0)
        self.context.fill()
      else:
        # ratio 0.4 looks better than 0.5 because the area difference
        # is percieved more than length difference
        ratio = 0.4
        xm1 = ratio*xa + (1-ratio)*x0 
        ym1 = ratio*ya + (1-ratio)*y0
        xm2 = (1-ratio)*(2*x2-x0) + ratio*(2*x1-xa)
        ym2 = (1-ratio)*(2*y2-y0) + ratio*(2*y1-ya)
        self.context.set_source_rgb( *color1)
        self._create_cairo_path( [(xa,ya), (xm1,ym1), (xm2,ym2), (2*x1-xa, 2*y1-ya)], closed=True)
        self.context.fill()
        self.context.set_source_rgb( *color2)
        self._create_cairo_path( [(xm1,ym1), (x0, y0), (2*x2-x0, 2*y2-y0), (xm2,ym2)], closed=True)
        self.context.fill()

    def draw_plain_or_colored_hatch( _start, _end):
      x1, y1 = _start
      x2, y2 = _end
      # no coloring now
      x, y, x0, y0 = geometry.find_parallel( x1, y1, x2, y2, self.wedge_width/2.0)
      xa, ya, xb, yb = geometry.find_parallel( x1, y1, x2, y2, self.line_width/2.0) 
      d = math.sqrt( (x1-x2)**2 + (y1-y2)**2) # length of the bond
      if d == 0:  
        return  # to prevent division by zero
      dx1 = (x0 - xa)/d 
      dy1 = (y0 - ya)/d 
      dx2 = (2*x2 -x0 -2*x1 +xa)/d 
      dy2 = (2*y2 -y0 -2*y1 +ya)/d 
      # we have to decide if the first line should be at the position of the first atom
      draw_start = 1  # is index not boolean
      if not v1 in self._vertex_to_bbox and v1.occupied_valency > 1:
        draw_start = 1
      draw_end = 1    # is added to index not boolean
      if not v2 in self._vertex_to_bbox and v2.occupied_valency > 1:
        draw_end = 0
      # adjust the step length
      step_size = 2*(self.line_width)
      ns = round( d / step_size) or 1
      step_size = d / ns
      # now we finally draw
      self.context.set_line_cap( cairo.LINE_CAP_BUTT)
      self.context.set_source_rgb( *color1)
      middle = 0.5 * (draw_start + int( round( d/ step_size)) + draw_end - 2)
      for i in range( draw_start, int( round( d/ step_size)) +draw_end):
        coords = [xa+dx1*i*step_size, ya+dy1*i*step_size, 2*x1-xa+dx2*i*step_size, 2*y1-ya+dy2*i*step_size] 
        if coords[0] == coords[2] and coords[1] == coords[3]:
          if (dx1+dx2) > (dy1+dy2): 
            coords[0] += 1
          else:
            coords[1] += 1
        self._create_cairo_path( [coords[:2],coords[2:]])
        if i >= middle:
          self.context.stroke()
          self.context.set_source_rgb( *color2)
      self.context.stroke()

    # code itself

    coords = self._where_to_draw_from_and_to( e)
    if not coords:
      return 
    start = coords[:2]
    end = coords[2:]
    v1, v2 = e.vertices
    color1 = self.atom_colors.get( v1.symbol, (0,0,0))
    color2 = self.atom_colors.get( v2.symbol, (0,0,0))
    has_shown_vertex = bool( [1 for _v in e.vertices if _v in self._vertex_to_bbox])

    if e.order == 1:
      if e.type == 'w':
        draw_plain_or_colored_wedge( start, end)
      elif e.type == 'h':
        draw_plain_or_colored_hatch( start, end)
      else:
        draw_plain_or_colored_line( start, end)
      
    if e.order == 2:
      side = 0
      # find how to center the bonds
      # rings have higher priority in setting the positioning
      in_ring = False
      for ring in self.molecule.get_smallest_independent_cycles():
        double_bonds = len( [b for b in self.molecule.vertex_subgraph_to_edge_subgraph(ring) if b.order == 2])
        if v1 in ring and v2 in ring:
          in_ring = True
          side += double_bonds * reduce( operator.add, [geometry.on_which_side_is_point( start+end, (a.x,a.y)) for a in ring if a!=v1 and a!=v2])
      # if rings did not decide, use the other neigbors
      if not side:
        for v in v1.neighbors + v2.neighbors:
          if v != v1 and v!= v2:
            side += geometry.on_which_side_is_point( start+end, (v.x, v.y))
      # if neighbors did not decide either
      if not side and in_ring:
        # we don't want centered bonds inside rings
        side = 1 # select arbitrary value
      if side:
        draw_plain_or_colored_line( start, end)
        x1, y1, x2, y2 = geometry.find_parallel( start[0], start[1], end[0], end[1], self.bond_width*misc.signum( side))
        # shorten the second line
        length = geometry.point_distance( x1,y1,x2,y2)
        if v2 not in self._vertex_to_bbox:
          x2, y2 = geometry.elongate_line( x1, y1, x2, y2, -self.bond_second_line_shortening*length)
        if v1 not in self._vertex_to_bbox:
          x1, y1 = geometry.elongate_line( x2, y2, x1, y1, -self.bond_second_line_shortening*length)
        draw_plain_or_colored_line( (x1, y1), (x2, y2), second=True)
      else:
        for i in (1,-1):
          x1, y1, x2, y2 = geometry.find_parallel( start[0], start[1], end[0], end[1], i*self.bond_width*0.5)
          draw_plain_or_colored_line( (x1, y1), (x2, y2))

    elif e.order == 3:
      self._draw_line( start, end, line_width=self.line_width)
      for i in (1,-1):
        x1, y1, x2, y2 = geometry.find_parallel( start[0], start[1], end[0], end[1], i*self.bond_width*0.7)
        draw_plain_or_colored_line( (x1, y1), (x2, y2), second=True)
    

  def _where_to_draw_from_and_to( self, b):
    # at first check if the bboxes are not overlapping
    atom1, atom2 = b.vertices
    x1, y1 = atom1.x, atom1.y
    x2, y2 = atom2.x, atom2.y
    bbox1 = self._vertex_to_bbox.get( atom1, None)
    bbox2 = self._vertex_to_bbox.get( atom2, None)
    if bbox1 and bbox2 and geometry.do_rectangles_intersect( bbox1, bbox2):
      return None
    # then we continue with computation
    if bbox1:
      x1, y1 = geometry.intersection_of_line_and_rect( (x1,y1,x2,y2), bbox1, round_edges=0)
    if bbox2:
      x2, y2 = geometry.intersection_of_line_and_rect( (x1,y1,x2,y2), bbox2, round_edges=0)
    if geometry.point_distance( x1, y1, x2, y2) <= 1.0:
      return None
    else:
      return (x1, y1, x2, y2)


  def _draw_vertex( self, v):
    if v.symbol != "C":
      x = v.x
      y = v.y
      text = v.symbol
      if self.show_hydrogens_on_hetero:
        if v.free_valency == 1:
          text += "H"
        elif v.free_valency > 1:
          text += "H<sub>%d</sub>" % v.free_valency

      if v.charge == 1:
        text += "<sup>+</sup>"
      elif v.charge == -1:
        text += "<sup>&#x2212;</sup>"
      elif v.charge > 1:
        text += "<sup>%d+</sup>" % v.charge
      elif v.charge < -1:
        text += "<sup>%d&#x2212;</sup>" % abs( v.charge)

      if self.color_atoms:
        color = self.atom_colors.get( v.symbol, (0,0,0))
      else:
        color = (0,0,0)
      bbox = self._draw_text( (x,y), text, center_first_letter=True, color=color)
      bbox = geometry.expand_rectangle( bbox, self.space_around_atom)
      self._vertex_to_bbox[v] = bbox


  ## ------------------------------ lowlevel drawing methods ------------------------------

  def _draw_colored_line( self, start, end, line_width=1, capstyle=cairo.LINE_CAP_BUTT,
                          start_color=(0,0,0), end_color=(0,0,0)):
    x1,y1 = start
    x2,y2 = end
    length = geometry.point_distance( x1,y1,x2,y2)
    xn2,yn2 = geometry.elongate_line( x1,y1,x2,y2, -0.5*length)
    line1 = [(x1,y1),(xn2,yn2)]
    xn1,yn1 = geometry.elongate_line( x2,y2,x1,y1, -0.5*length)
    line2 = [(xn1,yn1),(x2,y2)]
    self.context.set_line_cap( cairo.LINE_CAP_BUTT) # this is forced here
    self.context.set_line_width( line_width)
    for line,color in zip( [line1,line2], [start_color,end_color]):
      self.context.set_source_rgb( *color)
      self._create_cairo_path( line, closed=False)
      self.context.stroke()
    

  def _draw_line( self, start, end, line_width=1, capstyle=cairo.LINE_CAP_BUTT, color=(0,0,0)):
    self.context.set_source_rgb( *color)
    # cap style
    self.context.set_line_cap( capstyle)
    # line width
    self.context.set_line_width( line_width)
    # the path itself 
    cs = [start, end]
    self._create_cairo_path( cs, closed=False)
    # stroke it
    self.context.stroke()


  def _draw_text( self, xy, text, font_name=None, font_size=None, center_first_letter=False, color=(0,0,0)):
    import xml.sax
    from sets import Set
    class text_chunk:
      def __init__( self, text, attrs=None):
        self.text = text
        self.attrs = attrs or Set()

    class FtextHandler ( xml.sax.ContentHandler):
      def __init__( self):
        xml.sax.ContentHandler.__init__( self)
        self._above = []
        self.chunks = []
        self._text = ""
      def startElement( self, name, attrs):
        self._closeCurrentText()
        self._above.append( name)
      def endElement( self, name):
        self._closeCurrentText()
        self._above.pop( -1)
      def _closeCurrentText( self):
        if self._text:
          self.chunks.append( text_chunk( self._text, attrs = Set( self._above)))
          self._text = ""
      def characters( self, data):
        self._text += data

    # parse the text for markup
    handler = FtextHandler()
    try:
      xml.sax.parseString( "<x>%s</x>" % text, handler)
    except:
      chunks = [text_chunk( text)]
    else:
      chunks = handler.chunks

    if not font_name:
      font_name = self.font_name
    if not font_size:
      font_size = self.font_size

    # font properties
    self.context.select_font_face( font_name)
    self.context.set_font_size( font_size)
    self.context.set_line_width( 1.0)
    asc, desc, letter_height, _a, _b = self.context.font_extents()
    x, y = xy
    if center_first_letter:
      xbearing, ybearing, width, height, x_advance, y_advance = self.context.text_extents( chunks[0].text[0])
      x -= x_advance / 2.0
      y += height / 2.0
    
    self.context.new_path()
    x1 = round( x)
    bbox = None
    for chunk in chunks:
      y1 = y
      if "sup" in chunk.attrs:
        y1 -= asc / 2
        self.context.set_font_size( int( font_size * 0.8))
      elif "sub" in chunk.attrs:
        y1 += asc / 2
        self.context.set_font_size( int( font_size * 0.8))
      else:
        self.context.set_font_size( font_size)
      xbearing, ybearing, width, height, x_advance, y_advance = self.context.text_extents( chunk.text)
      # background
      if self.add_background_to_text:
        self.context.rectangle( x1+xbearing, y1+ybearing, width, height)
        self.context.set_source_rgb( *self.background_color)
        self.context.fill()
        #self.context.set_line_width( 3)
        #self.context.stroke()
      # text itself
      _bbox = [x1+xbearing, y1+ybearing, x1+xbearing+width, y1+ybearing+height]
      self._bboxes.append( _bbox)
      # store bbox for the first chunk only
      if not bbox:
        bbox = _bbox
      self.context.set_source_rgb( *color)
      self.context.move_to( x1, y1)
      self.context.show_text( chunk.text)
      #self.context.fill()
      x1 += x_advance
    return bbox



  def _draw_rectangle( self, coords, fill_color=(1,1,1)):
    #outline = self.paper.itemcget( item, 'outline')
    x1, y1, x2, y2 = coords
    self.context.set_line_join( cairo.LINE_JOIN_MITER)
    self.context.rectangle( x1, y1, x2-x1, y2-y1)
    self.context.set_source_rgb( *fill_color)
    self.context.fill_preserve()
    #self.context.set_line_width( width)
    #self.set_cairo_color( outline)
    self.context.stroke()

    
    
  def _create_cairo_path( self, points, closed=False):
    x, y = points.pop( 0)
    self.context.move_to( x, y)
    for (x,y) in points:
      self.context.line_to( x, y)
    if closed:
      self.context.close_path()

  def _get_bbox( self):
    bbox = list( self._bboxes[0])
    for _bbox in self._bboxes[1:]:
      x1,y1,x2,y2 = _bbox
      if x1 < bbox[0]:
        bbox[0] = x1
      if y1 < bbox[1]:
        bbox[1] = y1
      if x2 > bbox[2]:
        bbox[2] = x2
      if y2 > bbox[3]:
        bbox[3] = y2
    return bbox


def mol_to_png( mol, filename):
  c = cairo_out()
  c.mol_to_cairo( mol, filename)


if __name__ == "__main__":

##   import smiles

##   mol = smiles.text_to_mol( "c1cc(C#N)ccc1N", calc_coords=30)
##   mol_to_png( mol, "output.png")

  import inchi
  mol = inchi.text_to_mol( "1/C7H6O2/c8-7(9)6-4-2-1-3-5-6/h1-5H,(H,8,9)", include_hydrogens=False, calc_coords=30)
  mol_to_png( mol, "output.png")


