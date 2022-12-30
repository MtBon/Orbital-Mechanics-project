function g=derivata(n,f)
for i=1:n
    syms x;
    f=diff(f);
end
g=f;
end

