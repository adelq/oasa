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

import xml.dom.minidom as dom
import dom_extensions
import transform
import geometry
import math
import misc
import operator
import copy


class svg_out:

  show_hydrogens_on_hetero = False
  margin = 15
  line_width = 2
  bond_width = 6

  def __init__( self):
    pass

  def mol_to_svg( self, mol):
    self.document = dom.Document()
    self.top = dom_extensions.elementUnder( self.document,
                                            "svg",
                                            attributes=(("xmlns", "http://www.w3.org/2000/svg"),
                                                        ("version", "1.0")))
    self.top = dom_extensions.elementUnder( self.top, "g",
                                            attributes=(("stroke", "#000"),
                                                        ("stroke-width", "1.0")))


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

    self.molecule = mol

    for e in copy.copy( mol.edges):
      self._draw_edge( e)
    for v in mol.vertices:
      self._draw_vertex( v)

    return self.document


  def paper_to_canvas_coord( self, x):
    return x


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
      # rings have higher priority in setting the positioning
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
    

  def _draw_vertex( self, v):
    if v.symbol != "C":
      x = v.x - 5
      y = v.y + 6
      x1 = x
      x2 = x + 12
      y1 = y - 12
      y2 = y + 2
      text = v.symbol
      if self.show_hydrogens_on_hetero:
        if v.free_valency == 1:
          text += "H"
        elif v.free_valency > 1:
          text += "H%d" % v.free_valency

      if v.charge == 1:
        text += "+"
      elif v.charge == -1:
        text += "-"
      elif v.charge > 1:
        text += str( v.charge) + "+"
      elif v.charge < -1:
        text += str( v.charge)
        
      self._draw_rectangle( self.transformer.transform_4( (x1, y1, x2, y2)), fill_color="#fff")
      self._draw_text( self.transformer.transform_xy(x,y), text)


  def _draw_line( self, start, end, line_width=1, capstyle=""):
    x1, y1 = start
    x2, y2 = end
    line = dom_extensions.elementUnder(self.top, 'line',
                                        (( 'x1', str( x1)),
                                         ( 'y1', str( y1)),
                                         ( 'x2', str( x2)),
                                         ( 'y2', str( y2))))


  def _draw_text( self, xy, text, font_name="Arial", font_size=16):
    x, y = xy
    dom_extensions.textOnlyElementUnder( self.top, "text", text,
                                         (( "x", str( x)),
                                          ( "y", str( y)),
                                          ( "font-family", font_name),
                                          ( "font-size", str( font_size)),
                                          ( 'fill', "#000")))


  def _draw_rectangle( self, coords, fill_color="#fff", stroke_color="#fff"):
    x, y, x2, y2 = coords
    dom_extensions.elementUnder( self.top, 'rect',
                                 (( 'x', str( x)),
                                  ( 'y', str( y)),
                                  ( 'width', str( x2-x)),
                                  ( 'height', str( y2-y)),
                                  ( 'fill', fill_color),
                                  ( 'stroke', stroke_color)))
    

def mol_to_svg( mol, filename):
  c = svg_out()
  tree = c.mol_to_svg( mol)
  f = file( filename, "w")
  f.write( tree.toxml())
  f.close()


if __name__ == "__main__":

  import inchi
  mol = inchi.text_to_mol( "1/C7H6O2/c8-7(9)6-4-2-1-3-5-6/h1-5H,(H,8,9)", include_hydrogens=False, calc_coords=30)
  mol_to_svg( mol, "output.svg")


