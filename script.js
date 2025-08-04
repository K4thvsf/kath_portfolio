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

// Project Modal System
const projectData = {
    'phageseeker': {
        title: 'PhageSeeker',
        overview: 'Bioinformatics tool for detecting and analyzing bacteriophage sequences in bacterial genome assemblies. Uses BLAST database searches and pandas for data processing to identify viral sequences in FASTA format files.',
        methodology: `
            <h3>Implementation Approach</h3>
            <div class="modal-code">import pandas as pd
from Bio import SeqIO
# BLAST database search implementation
# Sequence alignment and analysis pipeline</div>
            
            <h3>Analysis Pipeline</h3>
            <p>1. <strong>Sequence Input Processing:</strong> Parse FASTA format bacterial genome assemblies and prepare sequences for analysis.</p>
            
            <p>2. <strong>BLAST Database Search:</strong> Perform sequence alignment against comprehensive bacteriophage database to identify potential viral sequences.</p>
            
            <p>3. <strong>Data Processing with pandas:</strong> Use pandas DataFrames for efficient manipulation and analysis of BLAST search results, including hit scoring and filtering.</p>
            
            <p>4. <strong>Results Analysis:</strong> Apply statistical thresholds and quality filters to identify high-confidence bacteriophage sequences within bacterial genomes.</p>
            
            <p>5. <strong>Report Generation:</strong> Generate comprehensive reports of identified viral sequences with detailed annotations and confidence scores.</p>
            
            <h3>Key Features</h3>
            <p>• Automated pipeline for bacteriophage detection in bacterial genomes</p>
            <p>• Integration with BLAST for sequence similarity searches</p>
            <p>• Efficient data handling using pandas for large genomic datasets</p>
            <p>• Comprehensive reporting system for identified viral sequences</p>
        `,
        technologies: ['Python', 'pandas', 'BLAST Database', 'Sequence Analysis', 'Data Processing', 'FASTA Analysis'],
        github: 'https://github.com/PhageSeeker/PHAGESEEKER'
    },
    
    'tcga-glioma': {
        title: 'TCGA Glioma Data Analysis',
        overview: 'Comprehensive statistical analysis of 1,122 brain lower grade glioma samples from TCGA using R and Bioconductor. Complete workflow from data cleaning through gene set analysis, including heatmap clustering, differential expression analysis, and IDH correlation studies.',
        methodology: `
            <h3>Data Preprocessing & Quality Control</h3>
            <div class="modal-code">library(omicsdata)
library(limma)
library(edgeR)
# Data cleaning and preprocessing
dataset <- fetch_tcga_dataset("LGG")
filtered_data <- filterByExpr(dataset)</div>
            
            <h3>Complete Analysis Workflow</h3>
            <p>1. <strong>Data Cleaning:</strong> Initial data quality assessment and removal of low-quality samples and genes.</p>
            
            <p>2. <strong>Exploratory Data Analysis:</strong> Evaluation of CG content bias to identify potential technical artifacts and batch effects.</p>
            
            <p>3. <strong>Filtering and Normalization:</strong> Applied TMM normalization and filtered low-expression genes to ensure robust downstream analysis.</p>
            
            <p>4. <strong>Inferential Statistics:</strong> Built design matrices and performed statistical modeling for differential expression analysis.</p>
            
            <p>5. <strong>Contrast Matrix & DGE Object:</strong> Constructed contrast matrices to define specific comparisons between glioma subtypes.</p>
            
            <p>6. <strong>Mean Gene Expression Plotting:</strong> Generated visualization of average expression levels across sample groups.</p>
            
            <p>7. <strong>Pairwise Comparisons:</strong> Performed statistical tests between all possible subtype combinations.</p>
            
            <p>8. <strong>FDR Correction:</strong> Applied false discovery rate correction to control for multiple testing bias.</p>
            
            <p>9. <strong>Differential Expression Analysis:</strong> Identified differentially expressed genes (DEGs) with statistical significance.</p>
            
            <p>10. <strong>P-value Distribution Analysis:</strong> Assessed the distribution of p-values to validate statistical model assumptions.</p>
            
            <p>11. <strong>Heatmap and Clustering:</strong> Generated hierarchical clustering heatmaps to visualize expression patterns and sample relationships.</p>
            
            <p>12. <strong>Gene Set Analysis:</strong> Comprehensive pathway analysis including:</p>
            <p>&nbsp;&nbsp;&nbsp;&nbsp;• Over-Representation Analysis (ORA) method</p>
            <p>&nbsp;&nbsp;&nbsp;&nbsp;• KEGG pathway enrichment analysis</p>
            <p>&nbsp;&nbsp;&nbsp;&nbsp;• REACTOME pathways analysis</p>
            <p>&nbsp;&nbsp;&nbsp;&nbsp;• Gene Ontology (GO) term enrichment</p>
            
            <p>13. <strong>IDH Status Correlation:</strong> Analyzed correlation between gene expression profiles and IDH mutation status, a key prognostic marker in gliomas.</p>
            
            <h3>Visualization & Results</h3>
            <p>• Comprehensive heatmaps showing gene expression clustering patterns</p>
            <p>• Statistical plots including p-value distributions and expression boxplots</p>
            <p>• Pathway enrichment visualizations for biological interpretation</p>
            <p>• IDH correlation analysis revealing molecular subtype characteristics</p>
            
            <h3>Project Heatmap Visualization</h3>
            <div style="text-align: center; margin: 1.5rem 0;">
                <img src="heatmap_glioma_project.png" alt="TCGA Glioma Expression Heatmap with Hierarchical Clustering" style="max-width: 100%; height: auto; border-radius: 8px; border: 1px solid var(--border-color); box-shadow: 0 4px 8px rgba(0,0,0,0.2);">
                <p style="font-size: 0.9rem; color: var(--text-light); margin-top: 0.5rem; font-style: italic;">
                    Hierarchical clustering heatmap showing gene expression patterns across glioma molecular subtypes (LGr1-4) with clear separation of sample groups based on expression profiles.
                </p>
            </div>
        `,
        technologies: ['R', 'Bioconductor', 'omicsdata', 'Statistical Analysis', 'Data Visualization', 'Survival Analysis', 'Exploratory Analysis', 'Data Cleaning']
    },
    
    'orf-prediction': {
        title: 'ORF Prediction Tool',
        overview: 'Developed a computational tool for predicting open reading frames (ORFs) in nucleotide sequences, demonstrated with BRCA2 gene analysis. The tool implements comprehensive start/stop codon detection with genetic code translation.',
        methodology: `
            <h3>Tool Implementation</h3>
            <div class="modal-code">def find_atg(seq):
    pos = seq.find("ATG")
    assert pos >= 0, "No ATG start codon found"
    return pos

def find_stop_codons(seq):
    stop_codons = ["UAA", "UAG", "UGA"]
    # Convert DNA to RNA
    rna_seq = seq.replace('T', 'U')
    # Find nearest stop codon
    return min_stop_position</div>
            
            <h3>Implementation Steps</h3>
            <p>1. <strong>Start codon detection:</strong> Locate ATG sequences within the input nucleotide sequence using string search algorithms with proper error handling for sequences without start codons.</p>
            
            <p>2. <strong>Stop codon identification:</strong> Implement detection for all three stop codons (UAA, UAG, UGA) with frame-aware searching to maintain proper reading frame alignment.</p>
            
            <p>3. <strong>Genetic code translation:</strong> Convert nucleotide triplets to amino acids using a comprehensive genetic code dictionary covering all 64 possible codons.</p>
            
            <p>4. <strong>ORF validation:</strong> Ensure predicted ORFs meet minimum length requirements and maintain proper reading frame throughout the sequence.</p>
            
            <h3>Example Analysis - BRCA2 Gene</h3>
            <div class="modal-code">seq = "AAAATGCCTATTGGATCCAAAGAGAGGCCAACATTTTTTGAAATTTTTAAGACACGCTGCAACAAAGCAGATTTAGGACCAATAAGTCTTTAAATGCAATAA"
pos = find_atg(seq)
orf_sequence = extract_orf(seq[pos:])
amino_acids = translate_to_protein(orf_sequence)</div>
        `,
        technologies: ['Python', 'Bioinformatics', 'Genetic Code', 'Sequence Analysis', 'BRCA2', 'String Processing', 'Tool Development']
    },
    
    'car-genetic-algorithm': {
        title: 'Genetic Algorithm Implementation',
        overview: 'Implementation of genetic algorithms for car design optimization developed as part of systems biology coursework. Features population-based evolution with selection, crossover, and mutation operations to optimize car performance across multiple generations.',
        methodology: `
            <h3>Genetic Algorithm Implementation</h3>
            <div class="modal-code">def initial_population(size):
    return [random_car() for _ in range(size)]

def generation(cars):
    scores = carsim.race(cars)
    top_cars = select_best(cars, scores, 0.3)
    new_generation = breed_population(top_cars)
    return new_generation</div>
            
            <h3>Car Genome Representation</h3>
            <p>Each car is encoded as a 15-element genome vector containing:</p>
            <p>• <strong>Frame length</strong> (1 element): Total car length between 2.0-6.0 meters</p>
            <p>• <strong>Upper/Lower profiles</strong> (10 elements): Height values at 5 equispaced points</p>
            <p>• <strong>Wheel specifications</strong> (4 elements): Position and radius for left/right wheels</p>
            
            <h3>Evolution Process</h3>
            <p>1. <strong>Population initialization:</strong> Generate 100 random car designs with constraint validation</p>
            
            <p>2. <strong>Fitness evaluation:</strong> Race cars on tracks and combine completion percentage with race time using the formula: fitness = fraction_covered × 1000 - elapsed_time</p>
            
            <p>3. <strong>Selection:</strong> Choose top 30% performers for breeding based on combined fitness score</p>
            
            <p>4. <strong>Crossover breeding:</strong> Generate offspring by combining genome segments from parent cars with random crossover points</p>
            
            <p>5. <strong>Mutation:</strong> Apply random mutations with 1% probability and 0.1 strength to maintain genetic diversity</p>
            
            <h3>Constraint Handling</h3>
            <div class="modal-code">def adjust_car(car):
    car[0] = np.clip(car[0], 2.0, 6.0)  # Frame length
    # Ensure total height > 0.5m and < 5.0m
    upper_shape, lower_shape = adjust_shapes(car[1], car[2])
    # Validate wheel positions and sizes
    return adjusted_car</div>
            
            <h3>Results After 30 Generations</h3>
            <p>• Successfully evolved cars capable of completing challenging off-road tracks</p>
            <p>• Demonstrated adaptation to different track profiles through track switching experiments</p>
            <p>• Achieved significant performance improvements through evolutionary optimization</p>
        `,
        technologies: ['Python', 'NumPy', 'Genetic Algorithms', 'Optimization', 'Data Analysis', 'Statistical Analysis', 'Algorithm Development']
    }
};

function openProjectModal(projectId) {
    const modal = document.getElementById('projectModal');
    const modalTitle = document.getElementById('modalTitle');
    const modalBody = document.getElementById('modalBody');
    
    const project = projectData[projectId];
    if (!project) return;
    
    modalTitle.textContent = project.title;
    
    modalBody.innerHTML = `
        <div class="modal-section">
            <h3>Overview</h3>
            <p>${project.overview}</p>
        </div>
        
        <div class="modal-section">
            <h3>Technologies Used</h3>
            <div class="modal-tech-grid">
                ${project.technologies.map(tech => 
                    `<span class="modal-tech-tag">${tech}</span>`
                ).join('')}
            </div>
        </div>
        
        <div class="modal-section">
            <h3>Methodology & Implementation</h3>
            ${project.methodology}
        </div>
        
        ${project.github ? `
        <div class="modal-section" style="text-align: center; margin-top: 2rem;">
            <a href="${project.github}" target="_blank" class="btn btn-primary" style="display: inline-flex; align-items: center; gap: 0.5rem;">
                <svg width="20" height="20" viewBox="0 0 24 24" fill="currentColor">
                    <path d="M12 0c-6.626 0-12 5.373-12 12 0 5.302 3.438 9.8 8.207 11.387.599.111.793-.261.793-.577v-2.234c-3.338.726-4.033-1.416-4.033-1.416-.546-1.387-1.333-1.756-1.333-1.756-1.089-.745.083-.729.083-.729 1.205.084 1.839 1.237 1.839 1.237 1.07 1.834 2.807 1.304 3.492.997.107-.775.418-1.305.762-1.604-2.665-.305-5.467-1.334-5.467-5.931 0-1.311.469-2.381 1.236-3.221-.124-.303-.535-1.524.117-3.176 0 0 1.008-.322 3.301 1.23.957-.266 1.983-.399 3.003-.404 1.02.005 2.047.138 3.006.404 2.291-1.552 3.297-1.23 3.297-1.30.653 1.653.242 2.874.118 3.176.77.84 1.235 1.911 1.235 3.221 0 4.609-2.807 5.624-5.479 5.921.43.372.823 1.102.823 2.222v3.293c0 .319.192.694.801.576 4.765-1.589 8.199-6.086 8.199-11.386 0-6.627-5.373-12-12-12z"/>
                </svg>
                View Complete Project on GitHub
            </a>
        </div>
        ` : ''}
    `;
    
    modal.style.display = 'block';
    document.body.style.overflow = 'hidden';
}

function closeProjectModal() {
    const modal = document.getElementById('projectModal');
    modal.style.display = 'none';
    document.body.style.overflow = 'auto';
}

// Close modal when clicking outside of it
window.onclick = function(event) {
    const modal = document.getElementById('projectModal');
    if (event.target === modal) {
        closeProjectModal();
    }
}

// Close modal with Escape key
document.addEventListener('keydown', function(event) {
    if (event.key === 'Escape') {
        closeProjectModal();
    }
});

// Export for potential external use
window.PortfolioApp = {
    DNAAnimationSystem,
    Navigation,
    ScrollAnimations,
    ContactForm,
    openProjectModal,
    closeProjectModal
};