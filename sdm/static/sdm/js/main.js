var app = new Vue({
    el: "#LRTstation",
    vuetify: new Vuetify(),
    data: () => ({
        step: 0,
        items: ['Substitution', 'Insertion', 'Deletion'],
        characterize: {
            sequence: null,
            mutation_type: null,
            mismatched_bases: null
        }
    }),
    methods: {
        prev() {
            this.step--;
        },
        next() {
            this.step++;
        },
        submit() {
            alert("Submitted");
        }
    }
});
