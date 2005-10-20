#--------------------------------------------------------------------------
#     This file is part of OASA - a free chemical python library
#     Copyright (C) 2003-2005 Beda Kosata <beda@zirael.org>

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
import transform
import geometry
import math
import misc
import operator



class cairo_out:


  _caps = {'butt': cairo.LINE_CAP_BUTT,
           'round': cairo.LINE_CAP_ROUND,
           'projecting': cairo.LINE_CAP_SQUARE}


  _font_remap = {'helvetica': 'Arial',
                 'times': 'Times New Roman'}

  _joins = {'round': cairo.LINE_JOIN_ROUND,
            'miter': cairo.LINE_JOIN_MITER,
            'bevel': cairo.LINE_JOIN_BEVEL}


  def __init__( self):
    self.margin = 15
    self.line_width = 2
    self.bond_width = 6

  def mol_to_cairo( self, mol, filename):
    x1, y1, x2, y2 = None, None, None, None
    for v in mol.vertices:
      if x1 == None or x1 > v.x:
        x1 = v.x
      if x2 == None or x2 < v.x:
        x2 = v.x
      if y1 == None or y1 > v.y:
        y1 = v.y
      if y2 == None or y2 < v.y:
        y2 = v.y

    w = int( x2 - x1 + 2*self.margin)
    h = int( y2 - y1 + 2*self.margin)

    self.transformer = transform.transform()
    self.transformer.set_move( -x1+self.margin, -y1+self.margin)

    self.surface = cairo.ImageSurface( cairo.FORMAT_ARGB32, w, h)

    self.context = cairo.Context( self.surface)
    self.context.set_source_rgb( 1, 1, 1)
    self.context.rectangle( 0, 0, w, h)
    self.context.fill()

    self.context.set_source_rgb( 0, 0, 0)

    self.molecule = mol

    for e in mol.edges:
      self._draw_edge( e)
    for v in mol.vertices:
      self._draw_vertext( v)

    self.context.show_page()
    self.surface.write_to_png( filename)
    self.surface.finish()



  def paper_to_canvas_coord( self, x):
    return x
    #dpi = self.paper.winfo_fpixels( '254m')/10.0
    #return dpi*x/72


  def prepare_dumb_transformer( self):
    tr = transform.transform()
    tr.set_scaling( self.paper_to_canvas_coord( 1))
    return tr


  def _draw_edge( self, e):
    v1, v2 = e.vertices
    start = self.transformer.transform_xy( v1.x, v1.y)
    end = self.transformer.transform_xy( v2.x, v2.y)

    if e.order == 1:
      self._draw_line( start, end, line_width=self.line_width)

    if e.order == 2:
      side = 0
      # find how to center the bonds
      # rings have higher prioriry in setting the positioning
      for ring in self.molecule.get_smallest_independent_cycles():
        if v1 in ring and v2 in ring:
          side += reduce( operator.add, [geometry.on_which_side_is_point( start+end, (self.transformer.transform_xy( a.x,a.y))) for a in ring if a!=v1 and a!=v2])
      # if rings did not decide, use the other neigbors
      if not side:
        for v in v1.neighbors + v2.neighbors:
          if v != v1 and v!= v2:
            side += geometry.on_which_side_is_point( start+end, (self.transformer.transform_xy( v.x, v.y)))
      if side:
        self._draw_line( start, end, line_width=self.line_width)
        x1, y1, x2, y2 = geometry.find_parallel( start[0], start[1], end[0], end[1], self.bond_width*misc.signum( side))
        self._draw_line( (x1, y1), (x2, y2), line_width=self.line_width)
      else:
        for i in (1,-1):
          x1, y1, x2, y2 = geometry.find_parallel( start[0], start[1], end[0], end[1], i*self.bond_width*0.5)
          self._draw_line( (x1, y1), (x2, y2), line_width=self.line_width)
        
        
    elif e.order == 3:
      self._draw_line( start, end, line_width=self.line_width)
      for i in (1,-1):
        x1, y1, x2, y2 = geometry.find_parallel( start[0], start[1], end[0], end[1], i*self.bond_width*0.7)
        self._draw_line( (x1, y1), (x2, y2), line_width=self.line_width)
    

  def _draw_vertext( self, v):
    if v.symbol != "C":
      x = v.x - 5
      y = v.y + 6
      x1 = x
      x2 = x + 12
      y1 = y - 12
      y2 = y + 2
      text = v.symbol
      if v.charge == 1:
        text += "+"
      elif v.charge == -1:
        text += "-"
      elif v.charge > 1:
        text += str( v.charge) + "+"
      elif v.charge < -1:
        text += str( v.charge)
        
      self._draw_rectangle( self.transformer.transform_4( (x1, y1, x2, y2)), fill_color=(1,1,1))
      self._draw_text( self.transformer.transform_xy(x,y), text)


  def _draw_line( self, start, end, line_width=1, capstyle=cairo.LINE_CAP_BUTT):
    #self.context.set_line_join( self._joins[ join])
    # color
    #self.set_cairo_color( self.paper.itemcget( item, 'fill'))
    # line width
    self.context.set_line_width( line_width)
    # the path itself 
    cs = [start, end]
    self._create_cairo_path( cs, closed=False)
    # stroke it
    self.context.stroke()



  def _draw_text( self, xy, text, font_name="Arial", font_size=16):
    # color
    #self.set_cairo_color( self.paper.itemcget( item, 'fill'))
    # helvetica which is often used does not work for me - therefor I use remap
    self.context.set_source_rgb( 0,0,0)

    self.context.select_font_face( font_name)

    self.context.set_font_size( font_size)
    #asc, desc, height, _a, _b = self.context.font_extents()
    #xbearing, ybearing, width, height, x_advance, y_advance = self.context.text_extents( text)

    x, y = xy

    self.context.new_path()
    self.context.move_to( x, y)
    #self.context.move_to( x1 - (width - x2 + x1)/2 - xbearing, y)
    self.context.text_path( text)
    self.context.fill()



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


