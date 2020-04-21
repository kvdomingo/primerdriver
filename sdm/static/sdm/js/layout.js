document.addEventListener('DOMContentLoaded', () => {
    var cl = new cloudinary.Cloudinary({
        cloud_name: 'kdphotography-assets',
        secure: true,
    });

    cl.responsive();
});
