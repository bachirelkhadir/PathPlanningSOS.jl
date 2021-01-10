from manim import *
from xml.dom import minidom
import string


def style_dict_from_string(style):
    if style in [None, ""]: return None
    elements = [s.split(":") for s in style.split(";")]
    elements = [e for e in elements if len(e) == 2]
    dict_elements = {key:value for key, value in elements}

    return dict_elements

def attribute_to_float(attr):
    stripped_attr = "".join([
        char for char in attr
        if char in string.digits + "." + "-"
    ])
    return float(stripped_attr)

def process_val_from_dict(key, D):
    if key not in D:
        return None
    v = D[key]
    
    if v is None:
        return None
    
    if type(v)==str: 
        if v.lower() == "none" or v == "":
            return None
        if v[0] == '#':
            return v.upper()
    
    return attribute_to_float(v)

def process_fill_stroke(element):
    style = style_dict_from_string(element.getAttribute("style"))
    
    if style is None: return None

    fill_color = process_val_from_dict("fill",style)
    opacity    = process_val_from_dict("fill-opacity",style)
    stroke_color = process_val_from_dict("stroke",style)
    stroke_width = process_val_from_dict("stroke-width",style)
    stroke_opacity = process_val_from_dict("stroke-opacity",style)
    
    if fill_color == "NONE" or fill_color=="": fill_color = None
    if stroke_color == "NONE" or stroke_color=="": stroke_color = None
    
    return fill_color, opacity, stroke_color, stroke_width, stroke_opacity

def extract_styles_from_elem(element):
    result = []
    if not isinstance(element, minidom.Element):
        return result
    if element.tagName in ['g', 'svg', 'symbol']:
        result += sum([extract_styles_from_elem(child) for child in element.childNodes],[])
    elif element.tagName in ['circle','rect','ellipse','path','polygon','polyline']:
        result.append(process_fill_stroke(element))
    return [r for r in result if r is not None]

def parse_styles(svg):
    doc = minidom.parse(svg.file_path)
    styles = []
    
    for svg_elem in doc.getElementsByTagName("svg"):
        styles += extract_styles_from_elem(svg_elem)
        
    doc.unlink()
    
    return styles

def color_svg_like_file(svgmobject):
    if not svgmobject.unpack_groups:
        raise Exception("Coloring groups not implemented yet!")
    styles = parse_styles(svgmobject)
    
    for i,(elem,style) in enumerate(zip(svgmobject,styles)):
        fc, alpha, sc, sw, salpha = style
        
        if alpha == 0. or fc is None:
            alpha = 0.
            fc = None
        
        if sw == 0. or sw is None or sc is None or salpha==0.:
            salpha = 0.
            sw = 0.
            sc = None
            
        svgmobject[i].set_fill(color=fc,opacity=alpha)
        svgmobject[i].set_stroke(color=sc,width=sw,opacity=salpha)
