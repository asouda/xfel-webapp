(() => {
  const canvas = document.getElementById("field");
  if (!canvas) return;
  const ctx = canvas.getContext("2d");
  const reduce = window.matchMedia("(prefers-reduced-motion: reduce)").matches;
  const darkColors = [
    [90, 150, 255],
    [70, 120, 255],
    [150, 90, 255],
    [190, 80, 255],
    [120, 160, 255],
  ];
  const lightColors = [
    [37, 84, 196],
    [46, 70, 190],
    [104, 58, 186],
    [132, 48, 176],
    [28, 96, 176],
  ];

  function colors() {
    return document.documentElement.getAttribute("data-theme") === "light" ? lightColors : darkColors;
  }
  let particles = [];

  function resize() {
    canvas.width = window.innerWidth;
    canvas.height = window.innerHeight;
    seed();
    draw();
  }

  function bell() {
    return (Math.random() + Math.random() + Math.random()) / 3;
  }

  function seed() {
    const count = Math.round(Math.min(110, Math.max(50, canvas.width * canvas.height / 18000)));
    particles = [];
    for (let i = 0; i < count; i += 1) {
      const large = Math.random() < 0.08;
      particles.push({
        x: bell() * canvas.width,
        y: bell() * canvas.height,
        r: large ? 18 + Math.random() * 46 : 1.2 + Math.random() * 3.4,
        vx: (Math.random() - 0.5) * (large ? 0.08 : 0.22),
        vy: (Math.random() - 0.5) * (large ? 0.08 : 0.22),
        a: large ? 0.05 + Math.random() * 0.08 : 0.25 + Math.random() * 0.55,
        colorIndex: i % colors().length,
      });
    }
  }

  function draw() {
    const palette = colors();
    const light = document.documentElement.getAttribute("data-theme") === "light";
    ctx.clearRect(0, 0, canvas.width, canvas.height);
    for (const p of particles) {
      const color = palette[p.colorIndex % palette.length];
      const alpha = light ? (p.r > 10 ? p.a * 1.8 : p.a * 0.7) : p.a;
      ctx.beginPath();
      ctx.fillStyle = `rgba(${color[0]}, ${color[1]}, ${color[2]}, ${alpha})`;
      ctx.shadowColor = `rgba(${color[0]}, ${color[1]}, ${color[2]}, ${Math.min(0.9, alpha + 0.2)})`;
      ctx.shadowBlur = p.r > 10 ? 40 : 16;
      ctx.arc(p.x, p.y, p.r, 0, Math.PI * 2);
      ctx.fill();
    }
    ctx.shadowBlur = 0;
  }

  function tick() {
    for (const p of particles) {
      p.x += p.vx;
      p.y += p.vy;
      if (p.x < -40) p.x = canvas.width + 40;
      if (p.x > canvas.width + 40) p.x = -40;
      if (p.y < -40) p.y = canvas.height + 40;
      if (p.y > canvas.height + 40) p.y = -40;
    }
    draw();
    requestAnimationFrame(tick);
  }

  window.addEventListener("xfel-theme", draw);
  window.addEventListener("resize", resize);
  resize();
  if (!reduce) requestAnimationFrame(tick);
})();
