import manim as m


def align_group_text(group, dir=m.LEFT):
    for g in group[1:]:
        g.align_to(group[0], dir)
    return group


def stack_group_text(group, dir=m.DOWN, buff=m.DEFAULT_MOBJECT_TO_MOBJECT_BUFFER):
    for g_prev, g in zip(group, group[1:]):
        g.next_to(g_prev, dir, buff=buff)
    return group
