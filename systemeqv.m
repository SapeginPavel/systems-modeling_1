function u = systemeqv(a, k)
    u = a + k * log(rand() / (1 - rand())); % генерация величины, распределенной по закону логистического распределения.
end
