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
    
    // Initialize language system with automatic detection
    new LanguageSystem();
    
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
        overview: {
            en: 'Bioinformatics tool for detecting and analyzing bacteriophage sequences in bacterial genome assemblies. Uses BLAST database searches and pandas for data processing to identify viral sequences in FASTA format files.',
            pt: 'Ferramenta de bioinformática para detectar e analisar sequências de bacteriófagos em montagens de genoma bacteriano. Utiliza buscas no banco BLAST e pandas para processamento de dados, identificando sequências virais em arquivos FASTA.'
        },
        methodology: {
            en: `
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
            pt: `
                <h3>Abordagem de Implementação</h3>
                <div class="modal-code">import pandas as pd
from Bio import SeqIO
# Implementação de busca no banco BLAST
# Pipeline de alinhamento e análise de sequências</div>
                
                <h3>Pipeline de Análise</h3>
                <p>1. <strong>Processamento de Sequências:</strong> Analisar montagens de genoma bacteriano em formato FASTA e preparar sequências para análise.</p>
                <p>2. <strong>Busca no Banco BLAST:</strong> Realizar alinhamento de sequências contra banco abrangente de bacteriófagos para identificar potenciais sequências virais.</p>
                <p>3. <strong>Processamento com pandas:</strong> Usar DataFrames pandas para manipulação e análise eficiente dos resultados BLAST, incluindo pontuação e filtragem.</p>
                <p>4. <strong>Análise de Resultados:</strong> Aplicar limites estatísticos e filtros de qualidade para identificar sequências de bacteriófagos de alta confiança.</p>
                <p>5. <strong>Geração de Relatórios:</strong> Gerar relatórios abrangentes de sequências virais identificadas com anotações detalhadas.</p>
                
                <h3>Características Principais</h3>
                <p>• Pipeline automatizado para detecção de bacteriófagos em genomas bacterianos</p>
                <p>• Integração com BLAST para buscas de similaridade de sequências</p>
                <p>• Manipulação eficiente de dados usando pandas para grandes datasets genômicos</p>
                <p>• Sistema abrangente de relatórios para sequências virais identificadas</p>
            `
        },
        technologies: ['Python', 'pandas', 'BLAST Database', 'Sequence Analysis', 'Data Processing', 'FASTA Analysis'],
        github: 'https://github.com/PhageSeeker/PHAGESEEKER'
    },
    
    'tcga-glioma': {
        title: 'TCGA Glioma Data Analysis',
        overview: {
            en: 'Comprehensive statistical analysis of 1,122 brain lower grade glioma samples from TCGA using R and Bioconductor. Complete workflow from data cleaning through gene set analysis, including heatmap clustering, differential expression analysis, and IDH correlation studies.',
            pt: 'Análise estatística completa de 1.122 amostras de glioma cerebral de baixo grau do TCGA usando R e Bioconductor. Workflow completo desde limpeza de dados até análise de conjuntos gênicos, incluindo clustering em heatmap, análise de expressão diferencial e estudos de correlação IDH.'
        },
        methodology: {
            en: `
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
            pt: `
                <h3>Pré-processamento & Controle de Qualidade</h3>
                <div class="modal-code">library(omicsdata)
library(limma)
library(edgeR)
# Limpeza e pré-processamento dos dados
dataset <- fetch_tcga_dataset("LGG")
filtered_data <- filterByExpr(dataset)</div>
                
                <h3>Workflow Completo de Análise</h3>
                <p>1. <strong>Limpeza de Dados:</strong> Avaliação inicial da qualidade e remoção de amostras e genes de baixa qualidade.</p>
                <p>2. <strong>Análise Exploratória:</strong> Avaliação de viés de conteúdo CG para identificar artefatos técnicos e efeitos de lote.</p>
                <p>3. <strong>Filtragem e Normalização:</strong> Aplicação de normalização TMM e filtragem de genes de baixa expressão.</p>
                <p>4. <strong>Estatística Inferencial:</strong> Construção de matrizes de design e modelagem estatística para análise de expressão diferencial.</p>
                <p>5. <strong>Matriz de Contraste:</strong> Construção de matrizes de contraste para comparações específicas entre subtipos de glioma.</p>
                <p>6. <strong>Gráficos de Expressão:</strong> Visualização dos níveis médios de expressão entre grupos de amostras.</p>
                <p>7. <strong>Comparações Pareadas:</strong> Testes estatísticos entre todas as combinações possíveis de subtipos.</p>
                <p>8. <strong>Correção FDR:</strong> Aplicação de correção de taxa de descoberta falsa para controlar viés de testes múltiplos.</p>
                <p>9. <strong>Análise de Expressão Diferencial:</strong> Identificação de genes diferencialmente expressos (DEGs) com significância estatística.</p>
                <p>10. <strong>Análise de Distribuição de P-valores:</strong> Avaliação da distribuição de p-valores para validar suposições do modelo estatístico.</p>
                <p>11. <strong>Heatmap e Clustering:</strong> Geração de heatmaps de clustering hierárquico para visualizar padrões de expressão.</p>
                <p>12. <strong>Análise de Conjuntos Gênicos:</strong> Análise abrangente de vias incluindo:</p>
                <p>&nbsp;&nbsp;&nbsp;&nbsp;• Método de Análise de Sobre-representação (ORA)</p>
                <p>&nbsp;&nbsp;&nbsp;&nbsp;• Análise de enriquecimento de vias KEGG</p>
                <p>&nbsp;&nbsp;&nbsp;&nbsp;• Análise de vias REACTOME</p>
                <p>&nbsp;&nbsp;&nbsp;&nbsp;• Enriquecimento de termos Gene Ontology (GO)</p>
                <p>13. <strong>Correlação Status IDH:</strong> Análise de correlação entre perfis de expressão gênica e status de mutação IDH.</p>
                
                <h3>Visualização & Resultados</h3>
                <p>• Heatmaps abrangentes mostrando padrões de clustering de expressão gênica</p>
                <p>• Gráficos estatísticos incluindo distribuições de p-valores e boxplots de expressão</p>
                <p>• Visualizações de enriquecimento de vias para interpretação biológica</p>
                <p>• Análise de correlação IDH revelando características de subtipos moleculares</p>
                
                <h3>Visualização do Heatmap do Projeto</h3>
                <div style="text-align: center; margin: 1.5rem 0;">
                    <img src="heatmap_glioma_project.png" alt="Heatmap de Expressão TCGA Glioma com Clustering Hierárquico" style="max-width: 100%; height: auto; border-radius: 8px; border: 1px solid var(--border-color); box-shadow: 0 4px 8px rgba(0,0,0,0.2);">
                    <p style="font-size: 0.9rem; color: var(--text-light); margin-top: 0.5rem; font-style: italic;">
                        Heatmap de clustering hierárquico mostrando padrões de expressão gênica entre subtipos moleculares de glioma (LGr1-4) com clara separação de grupos baseada em perfis de expressão.
                    </p>
                </div>
            `
        },
        technologies: ['R', 'Bioconductor', 'omicsdata', 'Statistical Analysis', 'Data Visualization', 'Survival Analysis', 'Exploratory Analysis', 'Data Cleaning']
    },
    
    'orf-prediction': {
        title: 'ORF Prediction Tool',
        overview: {
            en: 'Developed a computational tool for predicting open reading frames (ORFs) in nucleotide sequences, demonstrated with BRCA2 gene analysis. The tool implements comprehensive start/stop codon detection with genetic code translation.',
            pt: 'Desenvolvi uma ferramenta computacional para predizer quadros abertos de leitura (ORFs) em sequências de nucleotídeos, demonstrada com análise do gene BRCA2. A ferramenta implementa detecção abrangente de códons de início/parada com tradução do código genético.'
        },
        methodology: {
            en: `
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
            pt: `
                <h3>Implementação da Ferramenta</h3>
                <div class="modal-code">def find_atg(seq):
    pos = seq.find("ATG")
    assert pos >= 0, "Códon de início ATG não encontrado"
    return pos

def find_stop_codons(seq):
    stop_codons = ["UAA", "UAG", "UGA"]
    # Converter DNA para RNA
    rna_seq = seq.replace('T', 'U')
    # Encontrar códon de parada mais próximo
    return min_stop_position</div>
                
                <h3>Etapas de Implementação</h3>
                <p>1. <strong>Detecção de códon de início:</strong> Localizar sequências ATG dentro da sequência de nucleotídeos usando algoritmos de busca com tratamento de erros.</p>
                <p>2. <strong>Identificação de códons de parada:</strong> Implementar detecção para os três códons de parada (UAA, UAG, UGA) mantendo alinhamento correto do quadro de leitura.</p>
                <p>3. <strong>Tradução do código genético:</strong> Converter tripletos de nucleotídeos em aminoácidos usando dicionário completo do código genético.</p>
                <p>4. <strong>Validação de ORF:</strong> Garantir que ORFs preditos atendam requisitos mínimos de comprimento e mantenham quadro de leitura adequado.</p>
                
                <h3>Exemplo de Análise - Gene BRCA2</h3>
                <div class="modal-code">seq = "AAAATGCCTATTGGATCCAAAGAGAGGCCAACATTTTTTGAAATTTTTAAGACACGCTGCAACAAAGCAGATTTAGGACCAATAAGTCTTTAAATGCAATAA"
pos = find_atg(seq)
orf_sequence = extract_orf(seq[pos:])
amino_acids = translate_to_protein(orf_sequence)</div>
            `
        },
        technologies: ['Python', 'Bioinformatics', 'Genetic Code', 'Sequence Analysis', 'BRCA2', 'String Processing', 'Tool Development']
    },
    
    'car-genetic-algorithm': {
        title: 'Genetic Algorithm Implementation',
        overview: {
            en: 'Implementation of genetic algorithms for car design optimization developed as part of systems biology coursework. Features population-based evolution with selection, crossover, and mutation operations to optimize car performance across multiple generations.',
            pt: 'Implementação de algoritmos genéticos para otimização de design de carros desenvolvida como parte do curso de biologia de sistemas. Inclui evolução baseada em população com operações de seleção, cruzamento e mutação para otimizar performance de carros ao longo de múltiplas gerações.'
        },
        methodology: {
            en: `
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
                <p>2. <strong>Fitness evaluation:</strong> Race cars on tracks and combine completion percentage with race time</p>
                <p>3. <strong>Selection:</strong> Choose top 30% performers for breeding based on combined fitness score</p>
                <p>4. <strong>Crossover breeding:</strong> Generate offspring by combining genome segments from parent cars</p>
                <p>5. <strong>Mutation:</strong> Apply random mutations with 1% probability to maintain genetic diversity</p>
                
                <h3>Results After 30 Generations</h3>
                <p>• Successfully evolved cars capable of completing challenging off-road tracks</p>
                <p>• Demonstrated adaptation to different track profiles</p>
                <p>• Achieved significant performance improvements through evolutionary optimization</p>
            `,
            pt: `
                <h3>Implementação do Algoritmo Genético</h3>
                <div class="modal-code">def initial_population(size):
    return [random_car() for _ in range(size)]

def generation(cars):
    scores = carsim.race(cars)
    top_cars = select_best(cars, scores, 0.3)
    new_generation = breed_population(top_cars)
    return new_generation</div>
                
                <h3>Representação do Genoma do Carro</h3>
                <p>Cada carro é codificado como um vetor genoma de 15 elementos contendo:</p>
                <p>• <strong>Comprimento da estrutura</strong> (1 elemento): Comprimento total do carro entre 2,0-6,0 metros</p>
                <p>• <strong>Perfis superior/inferior</strong> (10 elementos): Valores de altura em 5 pontos equiespaçados</p>
                <p>• <strong>Especificações das rodas</strong> (4 elementos): Posição e raio para rodas esquerda/direita</p>
                
                <h3>Processo de Evolução</h3>
                <p>1. <strong>Inicialização da população:</strong> Gerar 100 designs aleatórios de carros com validação de restrições</p>
                <p>2. <strong>Avaliação de fitness:</strong> Correr carros em pistas e combinar porcentagem de conclusão com tempo de corrida</p>
                <p>3. <strong>Seleção:</strong> Escolher os 30% melhores performers para reprodução baseado no score de fitness</p>
                <p>4. <strong>Cruzamento reprodutivo:</strong> Gerar descendentes combinando segmentos genômicos de carros pais</p>
                <p>5. <strong>Mutação:</strong> Aplicar mutações aleatórias com probabilidade de 1% para manter diversidade genética</p>
                
                <h3>Resultados Após 30 Gerações</h3>
                <p>• Evolução bem-sucedida de carros capazes de completar pistas off-road desafiadoras</p>
                <p>• Demonstração de adaptação a diferentes perfis de pista</p>
                <p>• Melhorias significativas de performance através de otimização evolutiva</p>
            `
        },
        technologies: ['Python', 'NumPy', 'Genetic Algorithms', 'Optimization', 'Data Analysis', 'Statistical Analysis', 'Algorithm Development']
    }
};

function openProjectModal(projectId) {
    const modal = document.getElementById('projectModal');
    const modalTitle = document.getElementById('modalTitle');
    const modalBody = document.getElementById('modalBody');
    
    const project = projectData[projectId];
    if (!project) return;
    
    // Get current language from the LanguageSystem
    const currentLang = window.portfolioLanguage || 'en';
    
    modalTitle.textContent = project.title;
    
    // Get localized content
    const overview = project.overview[currentLang] || project.overview.en;
    const methodology = project.methodology[currentLang] || project.methodology.en;
    
    // Localized section titles
    const titles = {
        en: {
            overview: 'Overview',
            technologies: 'Technologies Used',
            methodology: 'Methodology & Implementation',
            github: 'View Complete Project on GitHub'
        },
        pt: {
            overview: 'Visão Geral',
            technologies: 'Tecnologias Utilizadas',
            methodology: 'Metodologia & Implementação',
            github: 'Ver Projeto Completo no GitHub'
        }
    };
    
    const t = titles[currentLang] || titles.en;
    
    modalBody.innerHTML = `
        <div class="modal-section">
            <h3>${t.overview}</h3>
            <p>${overview}</p>
        </div>
        
        <div class="modal-section">
            <h3>${t.technologies}</h3>
            <div class="modal-tech-grid">
                ${project.technologies.map(tech => 
                    `<span class="modal-tech-tag">${tech}</span>`
                ).join('')}
            </div>
        </div>
        
        <div class="modal-section">
            <h3>${t.methodology}</h3>
            ${methodology}
        </div>
        
        ${project.github ? `
        <div class="modal-section" style="text-align: center; margin-top: 2rem;">
            <a href="${project.github}" target="_blank" class="btn btn-primary" style="display: inline-flex; align-items: center; gap: 0.5rem;">
                <svg width="20" height="20" viewBox="0 0 24 24" fill="currentColor">
                    <path d="M12 0c-6.626 0-12 5.373-12 12 0 5.302 3.438 9.8 8.207 11.387.599.111.793-.261.793-.577v-2.234c-3.338.726-4.033-1.416-4.033-1.416-.546-1.387-1.333-1.756-1.333-1.756-1.089-.745.083-.729.083-.729 1.205.084 1.839 1.237 1.839 1.237 1.07 1.834 2.807 1.304 3.492.997.107-.775.418-1.305.762-1.604-2.665-.305-5.467-1.334-5.467-5.931 0-1.311.469-2.381 1.236-3.221-.124-.303-.535-1.524.117-3.176 0 0 1.008-.322 3.301 1.23.957-.266 1.983-.399 3.003-.404 1.02.005 2.047.138 3.006.404 2.291-1.552 3.297-1.23 3.297-1.30.653 1.653.242 2.874.118 3.176.77.84 1.235 1.911 1.235 3.221 0 4.609-2.807 5.624-5.479 5.921.43.372.823 1.102.823 2.222v3.293c0 .319.192.694.801.576 4.765-1.589 8.199-6.086 8.199-11.386 0-6.627-5.373-12-12-12z"/>
                </svg>
                ${t.github}
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

// Language Detection and Translation System
class LanguageSystem {
    constructor() {
        this.currentLanguage = 'en';
        this.translations = this.initTranslations();
        this.init();
    }
    
    initTranslations() {
        return {
            // Navigation
            'nav-home': { en: 'Home', pt: 'Início' },
            'nav-about': { en: 'About', pt: 'Sobre' },
            'nav-projects': { en: 'Projects', pt: 'Projetos' },
            'nav-contact': { en: 'Contact', pt: 'Contato' },
            
            // Hero Section
            'hero-greeting': { en: 'Hi, I\'m', pt: 'Oi, sou a' },
            'hero-title': { en: 'Bioinformatician & Data Scientist', pt: 'Bioinformata & Cientista de Dados' },
            'hero-description': { en: 'Passionate about transforming biological data into actionable insights through cutting-edge computational methods and machine learning techniques.', pt: 'Apaixonada por transformar dados biológicos complexos em descobertas usando computação avançada e machine learning.' },
            'hero-btn-work': { en: 'View My Work', pt: 'Ver Meus Projetos' },
            'hero-btn-contact': { en: 'Get In Touch', pt: 'Vamos Conversar' },
            
            // About Section
            'about-title': { en: 'About Me', pt: 'Sobre Mim' },
            'about-subtitle': { en: 'Bridging Biology and Technology', pt: 'Conectando Biologia e Tecnologia' },
            'about-p1': { en: 'I\'m a bioinformatician who loves turning complex biological data into discoveries that actually make a difference. I have an MSc in Molecular Biology with a focus on Bioinformatics and a BSc in Biomedicine. Currently living in Padua, Italy, I work on developing computational pipelines and machine learning models to better understand how biology works.', pt: 'Sou bioinformata e adoro transformar dados biológicos complexos em descobertas que realmente fazem diferença. Tenho mestrado em Biologia Molecular com foco em Bioinformática e graduação em Biomedicina. Moro em Pádua, na Itália, e trabalho desenvolvendo pipelines computacionais e modelos de machine learning para entender melhor como a biologia funciona.' },
            'about-p2': { en: 'I mainly work with Python, R, and SQL, always looking for creative ways to apply AI in biomedical research. One of the projects I\'m most proud of is PhageSeeker - a tool I co-created with colleagues to detect bacteriophages in bacterial genomes. I enjoy developing solutions that other researchers can actually use in their daily work.', pt: 'Trabalho principalmente com Python, R e SQL, sempre buscando formas criativas de aplicar IA na pesquisa biomédica. Um dos projetos que mais me orgulho é o PhageSeeker - uma ferramenta que criei junto com colegas para detectar bacteriófagos em genomas bacterianos. Gosto de desenvolver soluções que outros pesquisadores possam usar no dia a dia.' },
            'about-p3': { en: 'In my free time, you\'ll find me lost in a science fiction book, playing piano, planning my next trip, or curating the perfect indie playlists. I think these different passions help me think outside the box in research.', pt: 'Nos momentos livres, você vai me encontrar perdida em um livro de ficção científica, tocando piano, planejando a próxima viagem ou montando as melhores playlists indies. Acho que essas paixões diferentes me ajudam a pensar fora da caixa na pesquisa.' },
            
            // Skills
            'skills-programming': { en: 'Programming & Tools', pt: 'Programação & Ferramentas' },
            'skills-datascience': { en: 'Data Science & ML', pt: 'Ciência de Dados & ML' },
            'skills-bioinformatics': { en: 'Bioinformatics', pt: 'Bioinformática' },
            
            // Projects Section
            'projects-title': { en: 'Data Science Projects', pt: 'Projetos de Ciência de Dados' },
            'projects-subtitle': { en: 'Machine Learning & Computational Biology Solutions', pt: 'Soluções de Machine Learning & Biologia Computacional' },
            
            // Project 1 - PhageSeeker
            'project1-badge': { en: 'Bioinformatics • Data Analysis', pt: 'Bioinformática • Análise de Dados' },
            'project1-description': { en: 'Bioinformatics tool for detecting and analyzing bacteriophage sequences in bacterial genome assemblies. Uses BLAST database searches and pandas for data processing to identify viral sequences in FASTA format files.', pt: 'Ferramenta de bioinformática para detectar e analisar sequências de bacteriófagos em montagens de genoma bacteriano. Usa buscas no banco BLAST e pandas para processamento de dados para identificar sequências virais em arquivos formato FASTA.' },
            'project1-methods': { en: '<strong>Methods:</strong> BLAST database searches for sequence alignment, pandas for data manipulation and analysis, and comprehensive reporting of identified viral sequences in genomic datasets.', pt: '<strong>Métodos:</strong> Buscas no banco BLAST para alinhamento de sequências, pandas para manipulação e análise de dados, e relatórios abrangentes de sequências virais identificadas em datasets genômicos.' },
            
            // Project 2 - TCGA
            'project2-badge': { en: 'Data Analysis • Statistical Modeling', pt: 'Análise de Dados • Modelagem Estatística' },
            'project2-description': { en: 'Comprehensive statistical analysis of 1,122 brain lower grade glioma samples from TCGA using R and Bioconductor. Complete workflow from data cleaning through gene set analysis, including heatmap clustering, differential expression analysis, and IDH correlation studies.', pt: 'Análise estatística abrangente de 1.122 amostras de glioma cerebral de baixo grau do TCGA usando R e Bioconductor. Workflow completo desde limpeza de dados até análise de conjuntos de genes, incluindo clustering de heatmap, análise de expressão diferencial e estudos de correlação IDH.' },
            'project2-methods': { en: '<strong>Methods:</strong> RNA expression analysis, DNA methylation profiling, molecular subtype classification, and statistical analysis of clinical variables including age and ethnicity impacts on patient outcomes.', pt: '<strong>Métodos:</strong> Análise de expressão de RNA, perfil de metilação de DNA, classificação de subtipos moleculares e análise estatística de variáveis clínicas incluindo impactos de idade e etnia nos resultados dos pacientes.' },
            
            // Project 3 - ORF
            'project3-badge': { en: 'Bioinformatics • Data Science', pt: 'Bioinformática • Ciência de Dados' },
            'project3-description': { en: 'Python algorithm for predicting open reading frames (ORFs) in nucleotide sequences. Implements start codon detection, stop codon identification, and genetic code translation with BRCA2 gene analysis example.', pt: 'Algoritmo Python para predizer quadros de leitura aberta (ORFs) em sequências de nucleotídeos. Implementa detecção de códon de início, identificação de códon de parada e tradução de código genético com exemplo de análise do gene BRCA2.' },
            'project3-methods': { en: '<strong>Methods:</strong> ATG start codon detection, stop codon identification (UAA, UAG, UGA), and genetic code translation with error handling.', pt: '<strong>Métodos:</strong> Detecção de códon de início ATG, identificação de códons de parada (UAA, UAG, UGA) e tradução de código genético com tratamento de erros.' },
            
            // Project 4 - Genetic Algorithm
            'project4-badge': { en: 'Machine Learning • Optimization', pt: 'Machine Learning • Otimização' },
            'project4-description': { en: 'Machine learning implementation using genetic algorithms for optimization and predictive modeling. Features statistical analysis, fitness evaluation metrics, data-driven selection algorithms, and model performance optimization across 30+ generations.', pt: 'Implementação de machine learning usando algoritmos genéticos para otimização e modelagem preditiva. Apresenta análise estatística, métricas de avaliação de fitness, algoritmos de seleção baseados em dados e otimização de performance de modelo ao longo de mais de 30 gerações.' },
            'project4-methods': { en: '<strong>Methods:</strong> Population generation, fitness evaluation, selection algorithms, crossover breeding, and mutation operations for evolutionary optimization.', pt: '<strong>Métodos:</strong> Geração de população, avaliação de fitness, algoritmos de seleção, cruzamento reprodutivo e operações de mutação para otimização evolutiva.' },
            
            // Contact Section
            'contact-title': { en: 'Get In Touch', pt: 'Vamos Conversar' },
            'contact-subtitle': { en: 'Let\'s discuss how we can collaborate on exciting bioinformatics projects', pt: 'Que tal discutirmos como podemos colaborar em projetos incríveis de bioinformática?' },
            'contact-card-title': { en: 'Let\'s Connect', pt: 'Vamos nos Conectar' },
            'contact-card-text': { en: 'I\'m always interested in discussing new opportunities, research collaborations, or just chatting about the latest developments in bioinformatics and data science.', pt: 'Sempre tenho interesse em discutir novas oportunidades, colaborações de pesquisa ou simplesmente bater um papo sobre os últimos avanços em bioinformática e ciência de dados.' },
            'contact-email-title': { en: 'Contact Me', pt: 'Fale Comigo' },
            'contact-email-text': { en: 'Feel free to reach out to me directly via email:', pt: 'Sinta-se à vontade para entrar em contato diretamente por email:' },
            'contact-email-btn': { en: 'Send Email', pt: 'Enviar Email' },
            
            // Footer
            'footer-rights': { en: 'All rights reserved.', pt: 'Todos os direitos reservados.' },
            'footer-built': { en: 'Built with passion for bioinformatics.', pt: 'Feito com amor pela bioinformática.' }
        };
    }
    
    async init() {
        await this.detectUserLocation();
        this.createLanguageToggle();
        this.bindEvents();
        this.translatePage();
        
        // Set global variable for modals
        window.portfolioLanguage = this.currentLanguage;
    }
    
    async detectUserLocation() {
        try {
            // Try to get user's country from IP geolocation
            const response = await fetch('https://ipapi.co/json/');
            const data = await response.json();
            
            // Check if user is from a Portuguese-speaking country
            const portugueseCountries = ['BR', 'PT', 'AO', 'MZ', 'CV', 'GW', 'ST', 'TL', 'GQ'];
            
            if (portugueseCountries.includes(data.country_code)) {
                this.currentLanguage = 'pt';
            } else {
                this.currentLanguage = 'en';
            }
            
            console.log(`Language detected: ${this.currentLanguage} based on country: ${data.country_code}`);
        } catch (error) {
            console.log('Could not detect location, defaulting to English');
            this.currentLanguage = 'en';
        }
        
        // Check if user has previously selected a language
        const savedLanguage = localStorage.getItem('preferred-language');
        if (savedLanguage) {
            this.currentLanguage = savedLanguage;
        }
    }
    
    createLanguageToggle() {
        // Create language toggle if it doesn't exist
        let languageToggle = document.querySelector('.language-toggle');
        if (!languageToggle) {
            languageToggle = document.createElement('div');
            languageToggle.className = 'language-toggle';
            languageToggle.innerHTML = `
                <button id="lang-en" class="lang-btn ${this.currentLanguage === 'en' ? 'active' : ''}">EN</button>
                <button id="lang-pt" class="lang-btn ${this.currentLanguage === 'pt' ? 'active' : ''}">PT</button>
            `;
            document.body.appendChild(languageToggle);
        }
    }
    
    bindEvents() {
        const langEn = document.getElementById('lang-en');
        const langPt = document.getElementById('lang-pt');
        
        if (langEn) {
            langEn.addEventListener('click', () => this.setLanguage('en'));
        }
        if (langPt) {
            langPt.addEventListener('click', () => this.setLanguage('pt'));
        }
    }
    
    setLanguage(lang) {
        this.currentLanguage = lang;
        localStorage.setItem('preferred-language', lang);
        this.updateLanguageToggle();
        this.translatePage();
        
        // Update HTML lang attribute and global variable for modals
        document.documentElement.lang = lang;
        window.portfolioLanguage = lang;
    }
    
    updateLanguageToggle() {
        const langButtons = document.querySelectorAll('.lang-btn');
        langButtons.forEach(btn => {
            btn.classList.remove('active');
            if (btn.id === `lang-${this.currentLanguage}`) {
                btn.classList.add('active');
            }
        });
    }
    
    translatePage() {
        // Translate navigation
        this.translateElement('.nav-link[href="#home"]', 'nav-home');
        this.translateElement('.nav-link[href="#about"]', 'nav-about');
        this.translateElement('.nav-link[href="#projects"]', 'nav-projects');
        this.translateElement('.nav-link[href="#contact"]', 'nav-contact');
        
        // Translate hero section
        this.translateElement('.greeting', 'hero-greeting');
        this.translateElement('.hero-subtitle', 'hero-title');
        this.translateElement('.hero-description', 'hero-description');
        this.translateElement('.btn[href="#projects"]', 'hero-btn-work');
        this.translateElement('.btn[href="#contact"]', 'hero-btn-contact');
        
        // Translate about section
        this.translateElement('.about-section .section-title', 'about-title');
        this.translateElement('.about-section .section-subtitle', 'about-subtitle');
        
        // Translate about paragraphs
        const aboutParagraphs = document.querySelectorAll('.about-paragraph');
        if (aboutParagraphs.length >= 3) {
            aboutParagraphs[0].textContent = this.translations['about-p1'][this.currentLanguage];
            aboutParagraphs[1].textContent = this.translations['about-p2'][this.currentLanguage];
            aboutParagraphs[2].textContent = this.translations['about-p3'][this.currentLanguage];
        }
        
        // Translate skills categories
        const skillCategories = document.querySelectorAll('.skill-category h3');
        if (skillCategories.length >= 3) {
            skillCategories[0].textContent = this.translations['skills-programming'][this.currentLanguage];
            skillCategories[1].textContent = this.translations['skills-datascience'][this.currentLanguage];
            skillCategories[2].textContent = this.translations['skills-bioinformatics'][this.currentLanguage];
        }
        
        // Translate projects section
        this.translateElement('.projects-section .section-title', 'projects-title');
        this.translateElement('.projects-section .section-subtitle', 'projects-subtitle');
        
        // Translate individual projects
        this.translateProjectCards();
        
        // Translate contact section
        this.translateElement('.contact-section .section-title', 'contact-title');
        this.translateElement('.contact-section .section-subtitle', 'contact-subtitle');
        this.translateElement('.contact-card h3', 'contact-card-title');
        this.translateElement('.contact-card p', 'contact-card-text');
        this.translateElement('.contact-email-card h3', 'contact-email-title');
        this.translateElement('.contact-email-card p', 'contact-email-text');
        this.translateElement('.email-button', 'contact-email-btn');
        
        // Translate footer
        const footerParagraphs = document.querySelectorAll('.footer-text p');
        if (footerParagraphs.length >= 2) {
            footerParagraphs[0].innerHTML = `&copy; 2025 Kathleen Araujo. ${this.translations['footer-rights'][this.currentLanguage]}`;
            footerParagraphs[1].textContent = this.translations['footer-built'][this.currentLanguage];
        }
    }
    
    translateProjectCards() {
        const projectCards = document.querySelectorAll('.project-card');
        
        if (projectCards.length >= 4) {
            // Project 1 - PhageSeeker
            const badge1 = projectCards[0].querySelector('.data-science-badge');
            const desc1 = projectCards[0].querySelector('.project-description');
            const methods1 = projectCards[0].querySelector('.methods-text');
            
            if (badge1) badge1.textContent = this.translations['project1-badge'][this.currentLanguage];
            if (desc1) desc1.textContent = this.translations['project1-description'][this.currentLanguage];
            if (methods1) methods1.innerHTML = this.translations['project1-methods'][this.currentLanguage];
            
            // Project 2 - TCGA
            const badge2 = projectCards[1].querySelector('.data-science-badge');
            const desc2 = projectCards[1].querySelector('.project-description');
            const methods2 = projectCards[1].querySelector('.methods-text');
            
            if (badge2) badge2.textContent = this.translations['project2-badge'][this.currentLanguage];
            if (desc2) desc2.textContent = this.translations['project2-description'][this.currentLanguage];
            if (methods2) methods2.innerHTML = this.translations['project2-methods'][this.currentLanguage];
            
            // Project 3 - ORF
            const badge3 = projectCards[2].querySelector('.data-science-badge');
            const desc3 = projectCards[2].querySelector('.project-description');
            const methods3 = projectCards[2].querySelector('.methods-text');
            
            if (badge3) badge3.textContent = this.translations['project3-badge'][this.currentLanguage];
            if (desc3) desc3.textContent = this.translations['project3-description'][this.currentLanguage];
            if (methods3) methods3.innerHTML = this.translations['project3-methods'][this.currentLanguage];
            
            // Project 4 - Genetic Algorithm
            const badge4 = projectCards[3].querySelector('.data-science-badge');
            const desc4 = projectCards[3].querySelector('.project-description');
            const methods4 = projectCards[3].querySelector('.methods-text');
            
            if (badge4) badge4.textContent = this.translations['project4-badge'][this.currentLanguage];
            if (desc4) desc4.textContent = this.translations['project4-description'][this.currentLanguage];
            if (methods4) methods4.innerHTML = this.translations['project4-methods'][this.currentLanguage];
        }
    }
    
    translateElement(selector, translationKey) {
        const element = document.querySelector(selector);
        if (element && this.translations[translationKey]) {
            element.textContent = this.translations[translationKey][this.currentLanguage];
        }
    }
}

// Export for potential external use
window.PortfolioApp = {
    DNAAnimationSystem,
    Navigation,
    ScrollAnimations,
    ContactForm,
    LanguageSystem,
    openProjectModal,
    closeProjectModal
};