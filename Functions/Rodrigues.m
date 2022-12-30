function [v_rot] = Rodrigues(u,v,delta)

v_rot = v * cos(delta) + sin(delta) * cross(u,v) + dot(u,v) * u * (1 - cos(delta));