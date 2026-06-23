#pragma once

template <typename T, T v>
struct IntegralConstant
{
    static constexpr const T value = v;
    using value_type = T;
    using type = IntegralConstant;

    __device__ [[nodiscard]] inline consteval operator value_type() const noexcept
    {
        return value;
    }

    __device__ [[nodiscard]] inline consteval value_type operator()() const noexcept
    {
        return value;
    }
};

template <int Start, int End, typename F>
__device__ __forceinline__ constexpr void constexprFor(F &&f)
{
    if constexpr (Start < End)
    {
        f(IntegralConstant<int, Start>{});

        constexprFor<Start + 1, End>(
            static_cast<F &&>(f));
    }
}

template <
    typename Enum,
    int START,
    int END,
    class F>
__device__ void constexprEnumFor(F &&f)
{
    if constexpr (START < END)
    {
        f(IntegralConstant<
            Enum,
            static_cast<Enum>(START)>{});

        constexprEnumFor<
            Enum,
            START + 1,
            END>(f);
    }
}