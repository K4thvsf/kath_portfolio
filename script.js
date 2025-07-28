// Interactive Particles System - Design 2.1
class InteractiveParticles {
    constructor() {
        this.canvas = document.getElementById('interactive-particles');
        if (!this.canvas) {
            console.error('Interactive particles canvas not found!');
            return;
        }
        
        this.ctx = this.canvas.getContext('2d');
        this.particles = [];
        this.mouse = { x: window.innerWidth / 2, y: window.innerHeight / 2 };
        this.animationId = null;
        
        console.log('Interactive particles initializing...');
        this.init();
    }
    
    init() {
        this.setupCanvas();
        this.createParticles();
        this.bindEvents();
        this.animate();
    }
    
    setupCanvas() {
        this.resizeCanvas();
        window.addEventListener('resize', () => this.resizeCanvas());
    }
    
    resizeCanvas() {
        this.canvas.width = window.innerWidth;
        this.canvas.height = window.innerHeight;
    }
    
    createParticles() {
        const particleCount = Math.min(80, Math.floor((window.innerWidth * window.innerHeight) / 15000));
        
        for (let i = 0; i < particleCount; i++) {
            this.particles.push({
                x: Math.random() * this.canvas.width,
                y: Math.random() * this.canvas.height,
                baseX: Math.random() * this.canvas.width,
                baseY: Math.random() * this.canvas.height,
                vx: (Math.random() - 0.5) * 0.5,
                vy: (Math.random() - 0.5) * 0.5,
                size: Math.random() * 2 + 1,
                opacity: Math.random() * 0.3 + 0.1,
                color: '#67e8f9' // Design 2 accent color
            });
        }
    }
    
    bindEvents() {
        document.addEventListener('mousemove', (e) => {
            this.mouse.x = e.clientX;
            this.mouse.y = e.clientY;
        });
        
        window.addEventListener('resize', () => {
            this.particles = [];
            this.createParticles();
        });
    }
    
    animate() {
        this.ctx.clearRect(0, 0, this.canvas.width, this.canvas.height);
        
        this.particles.forEach((particle, index) => {
            // Calculate distance to mouse
            const dx = this.mouse.x - particle.x;
            const dy = this.mouse.y - particle.y;
            const distance = Math.sqrt(dx * dx + dy * dy);
            const maxDistance = 120;
            
            // Mouse interaction force
            if (distance < maxDistance) {
                const force = (maxDistance - distance) / maxDistance;
                const angle = Math.atan2(dy, dx);
                particle.x -= Math.cos(angle) * force * 2;
                particle.y -= Math.sin(angle) * force * 2;
            } else {
                // Return to base position when mouse is far away
                particle.x += (particle.baseX - particle.x) * 0.01;
                particle.y += (particle.baseY - particle.y) * 0.01;
            }
            
            // Add gentle floating movement
            particle.x += particle.vx;
            particle.y += particle.vy;
            
            // Boundary check and bounce
            if (particle.x <= 0 || particle.x >= this.canvas.width) {
                particle.vx *= -1;
                particle.baseX = particle.x;
            }
            if (particle.y <= 0 || particle.y >= this.canvas.height) {
                particle.vy *= -1;
                particle.baseY = particle.y;
            }
            
            // Draw particle
            this.ctx.beginPath();
            this.ctx.arc(particle.x, particle.y, particle.size, 0, Math.PI * 2);
            this.ctx.fillStyle = `rgba(103, 232, 249, ${particle.opacity})`;
            this.ctx.fill();
            
            // Draw connections to nearby particles
            this.particles.slice(index + 1).forEach(otherParticle => {
                const dx = particle.x - otherParticle.x;
                const dy = particle.y - otherParticle.y;
                const distance = Math.sqrt(dx * dx + dy * dy);
                
                if (distance < 100) {
                    const opacity = (100 - distance) / 100 * 0.1;
                    this.ctx.beginPath();
                    this.ctx.moveTo(particle.x, particle.y);
                    this.ctx.lineTo(otherParticle.x, otherParticle.y);
                    this.ctx.strokeStyle = `rgba(103, 232, 249, ${opacity})`;
                    this.ctx.lineWidth = 1;
                    this.ctx.stroke();
                }
            });
        });
        
        this.animationId = requestAnimationFrame(() => this.animate());
    }
    
    destroy() {
        if (this.animationId) {
            cancelAnimationFrame(this.animationId);
        }
    }
}

// Lightweight Blue DNA Helix Animation - Design 2 Optimized
class DNAAnimationSystem {
    constructor() {
        this.container = document.getElementById('dna-background');
        this.helixes = [];
        this.animationId = null;
        this.isReducedMotion = window.matchMedia('(prefers-reduced-motion: reduce)').matches;
        
        console.log('Lightweight DNA Animation System initializing...');
        
        if (!this.isReducedMotion && this.container) {
            this.init();
        }
    }
    
    init() {
        console.log('Lightweight DNA Animation System init called');
        this.createHelixes();
        this.bindEvents();
        this.animate();
    }
    
    createHelixes() {
        const isMobile = window.innerWidth <= 768;
        const helixCount = isMobile ? 8 : 12; // More helixes for better coverage
        
        console.log('Creating simple blue helixes:', { count: helixCount, isMobile });
        
        for (let i = 0; i < helixCount; i++) {
            this.createSimpleHelix(i);
        }
    }
    
    createSimpleHelix(index) {
        const helix = {
            element: document.createElement('div'),
            x: Math.random() * window.innerWidth,
            y: Math.random() * window.innerHeight,
            rotation: 0,
            rotationSpeed: 0.5 + Math.random() * 1.0,
            scale: 0.6 + Math.random() * 0.8,
            // Add movement properties
            moveX: (Math.random() - 0.5) * 0.5, // Slow horizontal drift
            moveY: (Math.random() - 0.5) * 0.3, // Slow vertical drift
            originalX: 0,
            originalY: 0
        };
        
        helix.element.className = 'simple-dna-helix';
        helix.element.innerHTML = `
            <svg width="80" height="200" viewBox="0 0 80 200" style="overflow: visible;">
                <defs>
                    <linearGradient id="helixGradient${index}" x1="0%" y1="0%" x2="100%" y2="100%">
                        <stop offset="0%" style="stop-color:#67e8f9;stop-opacity:0.8" />
                        <stop offset="50%" style="stop-color:#06b6d4;stop-opacity:1" />
                        <stop offset="100%" style="stop-color:#0891b2;stop-opacity:0.8" />
                    </linearGradient>
                    <filter id="glow${index}">
                        <feGaussianBlur stdDeviation="3" result="coloredBlur"/>
                        <feMerge> 
                            <feMergeNode in="coloredBlur"/>
                            <feMergeNode in="SourceGraphic"/>
                        </feMerge>
                    </filter>
                </defs>
                <!-- Simple helix path -->
                <path d="M20 0 Q40 25 20 50 Q0 75 20 100 Q40 125 20 150 Q0 175 20 200" 
                      stroke="url(#helixGradient${index})" 
                      stroke-width="3" 
                      fill="none" 
                      filter="url(#glow${index})"
                      opacity="0.7"/>
                <path d="M60 0 Q40 25 60 50 Q80 75 60 100 Q40 125 60 150 Q80 175 60 200" 
                      stroke="url(#helixGradient${index})" 
                      stroke-width="3" 
                      fill="none" 
                      filter="url(#glow${index})"
                      opacity="0.7"/>
                <!-- Connection lines -->
                <line x1="20" y1="25" x2="60" y2="25" stroke="#22d3ee" stroke-width="1.5" opacity="0.5"/>
                <line x1="20" y1="75" x2="60" y2="75" stroke="#22d3ee" stroke-width="1.5" opacity="0.5"/>
                <line x1="20" y1="125" x2="60" y2="125" stroke="#22d3ee" stroke-width="1.5" opacity="0.5"/>
                <line x1="20" y1="175" x2="60" y2="175" stroke="#22d3ee" stroke-width="1.5" opacity="0.5"/>
            </svg>
        `;
        
        // Store original positions for movement calculation
        helix.originalX = helix.x;
        helix.originalY = helix.y;
        
        helix.element.style.cssText = `
            position: absolute;
            left: ${helix.x}px;
            top: ${helix.y}px;
            transform: scale(${helix.scale});
            opacity: 0.4;
            pointer-events: none;
            z-index: 1;
        `;
        
        this.container.appendChild(helix.element);
        this.helixes.push(helix);
        
        console.log('Created simple helix:', index);
    }
    
    // Simple event binding - no complex interactions
    bindEvents() {
        window.addEventListener('resize', () => {
            this.handleResize();
        });
    }
    
    handleResize() {
        const isMobile = window.innerWidth <= 768;
        const targetHelixCount = isMobile ? 8 : 12;
        
        if (this.helixes.length !== targetHelixCount) {
            this.cleanup();
            this.createHelixes();
        }
    }
    
    // Simple animation loop - rotation and gentle movement
    animate() {
        this.helixes.forEach((helix, index) => {
            // Simple rotation on own axis - like in reference image
            helix.rotation += helix.rotationSpeed;
            
            // Add gentle floating movement
            helix.x += helix.moveX;
            helix.y += helix.moveY;
            
            // Keep helixes within screen bounds with gentle bounce
            if (helix.x > window.innerWidth + 100 || helix.x < -100) {
                helix.moveX *= -1;
            }
            if (helix.y > window.innerHeight + 100 || helix.y < -100) {
                helix.moveY *= -1;
            }
            
            // Apply rotation and position transforms
            helix.element.style.transform = `
                translate(${helix.x}px, ${helix.y}px)
                scale(${helix.scale})
                rotateZ(${helix.rotation}deg)
            `;
        });
        
        this.animationId = requestAnimationFrame(() => this.animate());
    }
    
    cleanup() {
        if (this.animationId) {
            cancelAnimationFrame(this.animationId);
        }
        this.helixes.forEach(helix => {
            if (helix.element.parentNode) {
                helix.element.parentNode.removeChild(helix.element);
            }
        });
        this.helixes = [];
    }
}

// Navigation functionality
class Navigation {
    constructor() {
        this.navbar = document.getElementById('navbar');
        this.navToggle = document.getElementById('nav-toggle');
        this.navMenu = document.getElementById('nav-menu');
        this.navLinks = document.querySelectorAll('.nav-link');
        
        this.init();
    }
    
    init() {
        this.bindEvents();
        this.updateActiveLink();
    }
    
    bindEvents() {
        // Mobile menu toggle
        this.navToggle.addEventListener('click', () => {
            this.navMenu.classList.toggle('active');
            this.navToggle.classList.toggle('active');
        });
        
        // Close mobile menu when clicking links
        this.navLinks.forEach(link => {
            link.addEventListener('click', () => {
                this.navMenu.classList.remove('active');
                this.navToggle.classList.remove('active');
            });
        });
        
        // Scroll spy for active navigation
        window.addEventListener('scroll', () => {
            this.handleScroll();
            this.updateActiveLink();
        });
    }
    
    handleScroll() {
        const scrolled = window.scrollY > 50;
        this.navbar.classList.toggle('scrolled', scrolled);
    }
    
    updateActiveLink() {
        const sections = document.querySelectorAll('section[id]');
        const scrollPos = window.scrollY + 100;
        
        sections.forEach(section => {
            const sectionTop = section.offsetTop;
            const sectionHeight = section.offsetHeight;
            const sectionId = section.getAttribute('id');
            
            if (scrollPos >= sectionTop && scrollPos < sectionTop + sectionHeight) {
                this.navLinks.forEach(link => {
                    link.classList.remove('active');
                    if (link.getAttribute('href') === `#${sectionId}`) {
                        link.classList.add('active');
                    }
                });
            }
        });
    }
}

// Scroll animations
class ScrollAnimations {
    constructor() {
        this.init();
    }
    
    init() {
        this.observeElements();
    }
    
    observeElements() {
        const observerOptions = {
            threshold: 0.1,
            rootMargin: '0px 0px -50px 0px'
        };
        
        const observer = new IntersectionObserver((entries) => {
            entries.forEach(entry => {
                if (entry.isIntersecting) {
                    entry.target.classList.add('visible');
                }
            });
        }, observerOptions);
        
        // Add fade-in class to elements that should animate
        const animateElements = document.querySelectorAll('.project-card, .skill-category, .contact-card');
        animateElements.forEach((el, index) => {
            el.classList.add('fade-in');
            el.style.transitionDelay = `${index * 0.1}s`;
            observer.observe(el);
        });
    }
}

// Contact functionality (no longer needed for form, but keeping class structure)
class ContactForm {
    constructor() {
        // No form to initialize anymore
        console.log('Contact section initialized - using mailto links');
    }
}

// Smooth scrolling for anchor links
function initSmoothScrolling() {
    document.querySelectorAll('a[href^="#"]').forEach(anchor => {
        anchor.addEventListener('click', function (e) {
            e.preventDefault();
            const target = document.querySelector(this.getAttribute('href'));
            if (target) {
                target.scrollIntoView({
                    behavior: 'smooth',
                    block: 'start'
                });
            }
        });
    });
}

// Performance optimization for mobile
function optimizeForMobile() {
    if (window.innerWidth <= 768) {
        // Reduce DNA animation complexity on mobile
        const dnaBackground = document.getElementById('dna-background');
        if (dnaBackground) {
            dnaBackground.style.opacity = '0.08';
        }
    }
}

// Preload critical resources
function preloadResources() {
    // Preload GitHub profile image
    const profileImg = new Image();
    profileImg.src = 'https://github.com/K4thvsf.png';
}

// Simple fallback DNA animation
function createSimpleDNAAnimation() {
    console.log('Using main DNA animation system - no fallback needed');
}

// Initialize everything when DOM is loaded
document.addEventListener('DOMContentLoaded', () => {
    console.log('DOM Content Loaded - Initializing portfolio');
    
    // Initialize core functionality
    new Navigation();
    new ScrollAnimations();
    new ContactForm();
    
    // Initialize DNA animation and interactive particles (only if motion is allowed)
    if (!window.matchMedia('(prefers-reduced-motion: reduce)').matches) {
        try {
            // Initialize DNA helixes
            new DNAAnimationSystem();
            
            // Initialize interactive particles
            new InteractiveParticles();
        } catch (error) {
            console.error('Animation systems failed, using fallback:', error);
            createSimpleDNAAnimation();
        }
    }
    
    // Initialize additional features
    initSmoothScrolling();
    optimizeForMobile();
    preloadResources();
    
    // Add loading animation completion
    document.body.classList.add('loaded');
});

// Handle page visibility changes to optimize performance
document.addEventListener('visibilitychange', () => {
    const dnaBackground = document.getElementById('dna-background');
    if (dnaBackground) {
        dnaBackground.style.animationPlayState = document.hidden ? 'paused' : 'running';
    }
});

// Export for potential external use
window.PortfolioApp = {
    DNAAnimationSystem,
    Navigation,
    ScrollAnimations,
    ContactForm
};